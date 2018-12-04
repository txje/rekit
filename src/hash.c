/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2017 Jeremy Wang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "klib/kvec.h" // C dynamic vector
#include "klib/khash.h" // C hash table/dictionary
#include "bnx.h"
#include "hash.h"
#include "dtw.h"
#include "chain.h"
#include "klib/ksort.h"

//#define aln_gt(a,b) ((a).score > (b).score)
//KSORT_INIT(aln_cmp, result, aln_gt)

uint8_t* get_fragments(label* labels, size_t n_labels, int bin_size, int rev) {
  uint8_t* frags = malloc(sizeof(uint8_t) * n_labels);
  if(n_labels == 0) {
    return frags;
  }
  frags[0] = labels[0].position / bin_size; // this *may* overflow if a fragment is >255 * bin_size, but it will just hash to val % 256
  int i;
  if(rev) {
    for(i = n_labels-1; i >= 1; i--) {
      frags[i] = (labels[i-1].position - labels[i].position) / bin_size;
    }
  } else {
    for(i = 1; i < n_labels; i++) {
      frags[i] = (labels[i].position - labels[i-1].position) / bin_size;
    }
  }
  return frags;
}

uint32_t* u32_get_fragments(label* labels, size_t n_labels, int bin_size, int rev) {
  uint32_t* frags = malloc(sizeof(uint32_t) * n_labels);
  if(n_labels == 0) {
    return frags;
  }
  frags[0] = labels[0].position / bin_size;
  int i;
  if(rev) {
    for(i = n_labels-1; i >= 1; i--) {
      frags[i] = (labels[i-1].position - labels[i].position) / bin_size;
    }
  } else {
    for(i = 1; i < n_labels; i++) {
      frags[i] = (labels[i].position - labels[i-1].position) / bin_size;
    }
  }
  return frags;
}


/*
 * nicks: an array of nick positions as produced by bntools (bn_file.c)
 * k: q-gram size
 * reverse: 0 - no, 1 - yes
 *
 * returns: 0 if successful, else 1
 */
int insert_rmap(label* labels, size_t n_labels, uint32_t read_id, int k, unsigned char reverse, khash_t(qgramHash) *db, int bin_size) {
  int i, j, l, absent, skip;

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  khint_t bin; // hash bin (result of kh_put)
  uint8_t* frags = get_fragments(labels, n_labels, bin_size, 0);
  if(n_labels < k) return 1;
  for(i = 0; i <= n_labels-k; i++) {
    khint_t qgram; // khint_t is probably u32
    //for(skip = 1; skip < k; skip++) {
      for(l = 0; l < pow(2, k-1); l++) { // iterate through a bit vector representing whether each position should be ceil'd
      //for(l = 0; l < 1; l++) { // iterate through a bit vector representing whether each position should be ceil'd
        //qgram = xratio_hash((labels+i), bin_size, skip);
        qgram = qgram_hash((frags+i), k, skip, l);
        //printf("# %dth %d-gram: %u\n", i, k, qgram);

        // insert qgram:readId,i into db
        bin = kh_put(qgramHash, db, qgram, &absent);
        if(absent) { // bin is empty (unset)
          kv_init(kh_value(db, bin));
        }
        // TODO: can reduce db memory usage by making this a pointer to a shared readPos among jittered hashes
        readPos r;
        r.readNum = (read_id << 1); // forward strand since its padded with a 0
        r.pos = i;
        kv_push(readPos, kh_value(db, bin), r);
      }
    //}
  }

  return 0;
}

// if reverse is true (1), forward and reverse signatures will be adjacent such that
// signature of read X forward is at index (2*X), and reverse is (2*X + 1)
void build_hash_db(cmap c, int k, khash_t(qgramHash) *db, int readLimit, int bin_size, int resolution_min) {

  label* filtered_labels;
  uint32_t f = 0;
  while (f < c.n_maps) {

    filtered_labels = malloc(c.molecules[f].n_labels * sizeof(label));
    int n_filtered_labels = filter_labels(c.molecules[f].labels, c.molecules[f].n_labels, filtered_labels, resolution_min);
    //int res = insert_rmap(c.labels[f], c.map_lengths[f], f, k, 0, db, bin_size); // forward strand only right now
    int res = insert_rmap(filtered_labels, n_filtered_labels, f, k, 0, db, bin_size); // forward strand only right now
    free(filtered_labels);
    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
}

// floor: the index of the fragment to floor
int jitter_bins(uint8_t *frags, int i, int k, int skip, khash_t(qgramHash) *db, khash_t(matchHash) *hits, int max_qgrams) {
  int j, m, absent;
  uint32_t l;
  khint_t bin; // hash bin (result of kh_put)

  //qgram = xratio_hash((labels+i), bin_size, skip, floor);
  //for(l = 0; l < pow(2, k-1); l++) { // iterate through a bit vector representing whether each position should be ceil'd
  for(l = 0; l < 1; l++) { // iterate through a bit vector representing whether each position should be ceil'd
    khint_t qgram = qgram_hash((frags+i), k, skip, l);

    //fprintf(stderr, "# %dth %d-gram (%u): %u\n", i, k, l, qgram);

    khint_t size_close, close;
    close = qgram;
    /*
    for(size_close = qgram > bin_size ? qgram-bin_size : qgram; size_close < qgram + bin_size + 1; size_close += bin_size) {
      for(close = size_close > 0 ? size_close-1 : 0; close < size_close+2; close++) {
    */
        //fprintf(stderr, "adding hash val %u\n", close);
        bin = kh_get(qgramHash, db, close);
        if(bin == kh_end(db)) // key not found, IDK what happens if you don't test this
          continue;
        matchVec matches = kh_val(db, bin);
        if(kv_size(matches) > max_qgrams) { // repetitive, ignore it
          continue;
        }
        for(m = 0; m < kv_size(matches); m++) {
          bin = kh_put(matchHash, hits, kv_A(matches, m).readNum>>1, &absent); // >>1 removes the fw/rv bit, which is always fw(0) right now
          if(absent) { // bin is empty (unset)
            kv_init(kh_value(hits, bin));
          }
          // check back to be sure we haven't already found this pair (this is possible now with the 4/5-mers)
          char exists = 0;
          for(j = kv_size(kh_value(hits, bin))-1; j > 0 && kv_A(kh_value(hits, bin), j).qpos == i; j--) {
            if(kv_A(kh_value(hits, bin), j).tpos == kv_A(matches, m).pos) {
              exists = 1;
              break;
            }
          }
          if(!exists) { // this pair was not already found
            posPair ppair; // to store matching query/target positions
            ppair.qpos = i;
            ppair.tpos = kv_A(matches, m).pos;
            kv_push(posPair, kh_value(hits, bin), ppair);
          }
        }
    /*
      }
    }
    */
  } // </for l in range(2)>
  return 0;
}

khash_t(matchHash)* lookup(label* labels, size_t n_labels, uint32_t read_id, int k, uint8_t rev, khash_t(qgramHash) *db, int max_qgrams, int bin_size) {
  int i;

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  khash_t(matchHash) *hits = kh_init(matchHash); // allocate hash table
  uint8_t* frags = get_fragments(labels, n_labels, bin_size, rev);
  if(n_labels < k) return hits;
  for(i = 0; i <= n_labels-k; i++) {
    int skip;
    //for(skip = 1; skip < k; skip++) {
    //for(skip = 1; skip <= 1; skip++) {
      int res = jitter_bins(frags, i, k, skip, db, hits, max_qgrams);
    //}
  }

  return hits;
}

void query_db(cmap b, int k, khash_t(qgramHash) *db, cmap c, FILE* o, int readLimit, int max_qgrams, int chain_threshold, float dtw_threshold, int bin_size, int min_labels, int start_mol, int end_mol) {
  int i, j, l, a;
  uint32_t target;
  int max_chains = 10000000; // this can be a parameter
  int max_alignments = 3; // this can be a parameter
  int match_score = 4; // ? idk what to do with this
  int max_gap = 50; // need to test/refine this
  int min_chain_length = 3; // need to test/refine this

  uint32_t f = start_mol;
  //label* filtered_labels;
  while (f <= end_mol) {
    if(b.molecules[f].n_labels < min_labels) {f++; continue;}; // enforce minimum number of labels to attempt alignment

    uint8_t qrev;

    kvec_t(result) alignments;
    kv_init(alignments);
    //result* alignments = malloc(max_chains * 2 * sizeof(result)); // we can actually get n_chains from each direction (fw/rv)
    //a = 0;

    for(qrev = 0; qrev <= 1; qrev++) {

      // ------ we do this kind of filtering in the simulation now ------
      //filtered_labels = malloc(b.map_lengths[f] * sizeof(label));
      //int n_filtered_labels = filter_labels(b.labels[f], b.map_lengths[f], filtered_labels, 500);
      fprintf(stderr, "# Hashing fragment of size %d with %d nicks\n", b.molecules[f].length, b.molecules[f].n_labels);
      khash_t(matchHash) *hits = lookup(b.molecules[f].labels, b.molecules[f].n_labels, f, k, qrev, db, max_qgrams, bin_size); // forward strand only right now
      //khash_t(matchHash) *hits = lookup(filtered_labels, n_filtered_labels, f, k, qrev, db, max_qgrams, bin_size); // forward strand only right now
      //free(filtered_labels);

      chain* chains = do_chain(hits, max_chains, match_score, max_gap, min_chain_length);
      int n_chains = 0;
      // count chains (if fewer than max_chains, the chains array is terminated by an empty chain vector)
      for(n_chains = 0; n_chains < max_chains && kv_size(chains[n_chains].anchors) > 0; n_chains++);

      int* starts = malloc(n_chains * sizeof(int));
      int* ends = malloc(n_chains * sizeof(int));
      uint32_t* refs = malloc(n_chains * sizeof(uint32_t));

      /*
      fprintf(stderr, "%d chains found with anchor sizes: ", n_chains);
      for(i = 0; i < n_chains; i++) {
        fprintf(stderr, "%d, ", kv_size(chains[i].anchors));
      }
      fprintf(stderr, "\n");
      */

      l = 0; // count of non-overlapping chains
      int last; // index of the last range that was merged
      for(j = 0; j < n_chains; j++) {
        if(kv_size(chains[j].anchors) < chain_threshold) continue;
        target = chains[j].ref; // the target is encoded in the chained score struct, do_chain() should have enforced that all chained anchors are from the same target

        // dynamic time warping
        // extract ref labels - expand bounds to encompass unmatched labels within query range
        int rst = kv_A(chains[j].anchors, 0).tpos;
        // esimated start position on ref is (anchor[0]_ref_pos - anchor[0]_query_pos)
        int est_rst = c.molecules[target].labels[rst].position - b.molecules[f].labels[kv_A(chains[j].anchors, 0).qpos].position;
        while(rst > 0 && c.molecules[target].labels[rst].position > est_rst)
          rst--;
        int ren = kv_A(chains[j].anchors, kv_size(chains[j].anchors)-1).tpos;
        // estimated end position on ref is (anchor[n]_ref_pos + (query_length - anchor[n]_query_pos))
        int est_ren = c.molecules[target].labels[kv_A(chains[j].anchors, kv_size(chains[j].anchors)-1).tpos].position + (b.molecules[f].labels[b.molecules[f].n_labels-1].position - b.molecules[f].labels[kv_A(chains[j].anchors, kv_size(chains[j].anchors)-1).qpos].position);
        while(ren < c.molecules[target].n_labels-1 && c.molecules[target].labels[ren].position < est_ren)
          ren++;

        /*
        fprintf(stderr, "chain %d\n", j);
        fprintf(stderr, "score %d\n", chains[j].score);
        fprintf(stderr, "est ref pos %d - %d\n", est_rst, est_ren);
        fprintf(stderr, "r indices %d - %d\n", rst, ren);
        */

        // loop through previous chain bounds and merge if they overlap
        last = -1;
        for(i = 0; i < l; i++) {
          if(target == refs[i] && ((last > -1 && starts[last] <= ends[i] && ends[last] >= starts[i]) || (last == -1 && rst <= ends[i] && ren >= starts[i]))) {
            if(last > -1) {
              refs[last] = -1; // unset this one since it was merged down
              starts[i] = starts[last] < starts[i] ? starts[last] : starts[i];
              ends[i] = ends[last] > ends[i] ? ends[last] : ends[i];
            } else {
              starts[i] = rst < starts[i] ? rst : starts[i];
              ends[i] = ren > ends[i] ? ren : ends[i];
            }
            // we can't stop here, we have to keep merging down
            last = i;
            //fprintf(stderr, "overlaps %d: %d-%d\n", refs[i], starts[i], ends[i]);
          }
        }
        // didn't overlap any
        if(last == -1) {
          starts[i] = rst;
          ends[i] = ren;
          refs[i] = target;
          l++;
          //fprintf(stderr, "new %d: %d-%d\n", refs[i], starts[i], ends[i]);
        }
      }
      n_chains = l; // includes those that were merged overlaps (ref == -1)

      for(j = 0; j < n_chains; j++) {
        if(refs[j] == -1) continue; // merged down
        // get fragment distances for DTW (no discretization)
        uint32_t* qfrags = u32_get_fragments(b.molecules[f].labels, b.molecules[f].n_labels, 1, 0); // get the fw ordered fragments, the reversal will be handled by the DTW
        uint32_t* rfrags = u32_get_fragments(c.molecules[refs[j]].labels+starts[j], ends[j]-starts[j]+1, 1, 0);
        //fprintf(stderr, "running dtw for read %d to ref %u %u-%u (of %u)\n", f, refs[j], starts[j], ends[j], c.map_lengths[refs[j]]);
        result aln = dtw(qfrags, rfrags, b.molecules[f].n_labels, ends[j]-starts[j]+1, -1, -1, 0.2, qrev); // ins_score, del_score, neutral_deviation
        aln.tstart += starts[j];
        aln.tend += starts[j];
        aln.ref = refs[j];
        //alignments[a++] = aln;
        kv_push(result, alignments, aln);
        free(qfrags);
        free(rfrags);
        if(aln.failed) {
          //fprintf(stderr, "q %d : ref %d DTW failed -- this should never happen\n", f, refs[j]);
          aln.score = -1; // to make sure it's sorted to the bottom
          continue;
        }
      }

      free(starts);
      free(ends);
      free(refs);
    } // </qrev>


    // sort alignments by (DTW) score decreasing
    ks_mergesort(aln_cmp, kv_size(alignments), alignments.a, 0);

    for(j = 0; j < (max_alignments < kv_size(alignments) ? max_alignments : kv_size(alignments)); j++) {
      if(kv_A(alignments, j).failed || kv_A(alignments, j).score < dtw_threshold) {
        fprintf(o, "%u\t-\t-\t-\t-\t%u\t-\t-\t%u\t-\t-\t-\t-\t-\t-\t-\t-\n", b.molecules[f].id, b.molecules[f].n_labels, b.molecules[f].length);
        continue;
      }

      // print chain output only
      /*
      printf("%d,%d,%d,%d", f, qrev, target, kv_size(chains[j].anchors));
      for(i = 0; i < kv_size(chains[j].anchors); i++)
        printf(",%u(%u):%u(%u)", kv_A(chains[j].anchors, i).qpos, b.labels[f][kv_A(chains[j].anchors, i).qpos].position, kv_A(chains[j].anchors, i).tpos, c.labels[target][kv_A(chains[j].anchors, i).tpos].position);
      printf("\n");
      */

      fprintf(o, "%u\t", b.molecules[f].id); // query id
      fprintf(o, "%u\t", c.molecules[kv_A(alignments, j).ref].id); // target id
      fprintf(o, "%u\t", qrev); // query reverse?
      fprintf(o, "%u\t", kv_A(alignments, j).qstart); // query start idx
      fprintf(o, "%u\t", kv_A(alignments, j).qend); // query end idx
      fprintf(o, "%u\t", b.molecules[f].n_labels); // query len idx
      fprintf(o, "%u\t", b.molecules[f].labels[kv_A(alignments, j).qstart].position); // query start
      fprintf(o, "%u\t", b.molecules[f].labels[kv_A(alignments, j).qend-1].position); // query end
      fprintf(o, "%u\t", b.molecules[f].length); // query len
      fprintf(o, "%u\t", kv_A(alignments, j).tstart); // ref start idx
      fprintf(o, "%u\t", kv_A(alignments, j).tend); // ref end idx
      fprintf(o, "%u\t", c.molecules[kv_A(alignments, j).ref].n_labels); // ref len idx
      fprintf(o, "%u\t", c.molecules[kv_A(alignments, j).ref].labels[kv_A(alignments, j).tstart].position); // ref start
      fprintf(o, "%u\t", c.molecules[kv_A(alignments, j).ref].labels[kv_A(alignments, j).tend-1].position); // ref end
      fprintf(o, "%u\t", c.molecules[kv_A(alignments, j).ref].length); // ref len
      fprintf(o, "%f\t", kv_A(alignments, j).score); // dtw score
      // dtw path
      for(i = 0; i < kv_size(kv_A(alignments, j).path); i++) {
        fprintf(o, "%c", kv_A(kv_A(alignments, j).path, i) == 0 ? '.' : (kv_A(kv_A(alignments, j).path, i) == 1 ? 'I' : 'D'));
      }
      fprintf(o, "\n");
    }

    if(readLimit > 0 && f >= readLimit) {
      break;
    }

    f++;
  }
}


/*
 * readLimit: maximum reads to process for BOTH database and query
 */
int hash_cmap(cmap b, cmap c, FILE* o, int q, int chain_threshold, float dtw_threshold, int max_qgrams, int readLimit, int bin_size, int resolution_min, int min_labels, int start_mol, int end_mol) {

  // ------------------------- Create hash database -----------------------------

  time_t t0 = time(NULL);

  // read BNX file, construct hash, including only forward direction
  // -- maybe assess doing this the opposite way at a later date, I'm not sure which will be faster
  khash_t(qgramHash) *db = kh_init(qgramHash);

  fprintf(stderr, "# Hashing %d cmap fragments\n", c.n_maps);
  build_hash_db(c, q, db, readLimit, bin_size, resolution_min);

  time_t t1 = time(NULL);
  fprintf(stderr, "# Hashed rmaps in %d seconds\n", (t1-t0));
  t0 = t1;
  // -------------------------------------------------------------------------------

  // ---------------------------- Look up queries in db ------------------------------
  fprintf(stderr, "# Querying %d bnx fragments\n", b.n_maps);
  query_db(b, q, db, c, o, readLimit, max_qgrams, chain_threshold, dtw_threshold, bin_size, min_labels, start_mol, end_mol);

  t1 = time(NULL);
  fprintf(stderr, "# Queried and output in %d seconds\n", (t1-t0));
  // ----------------------------------------------------------------------------------------

  // TODO: clean up vectors in each hash bin
  kh_destroy(qgramHash, db);
  return 0;
}
