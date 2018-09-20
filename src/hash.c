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
#include "chain.c"
#include "klib/ksort.h"

#define pos_pair_lt(a,b) ((a).tpos < (b).tpos)
KSORT_INIT(pos_pair_cmp, posPair, pos_pair_lt)

uint8_t* get_fragments(label* labels, size_t n_labels, int bin_size) {
  uint8_t* frags = malloc(sizeof(uint8_t) * n_labels);
  if(n_labels == 0) {
    return frags;
  }
  frags[0] = labels[0].position / bin_size; // this *may* overflow if a fragment is >255 * bin_size, but it will just hash to val % 256
  int i;
  for(i = 1; i < n_labels; i++) {
    frags[i] = (labels[i].position - labels[i-1].position) / bin_size;
  }
  return frags;
}

uint32_t* u32_get_fragments(label* labels, size_t n_labels, int bin_size) {
  uint32_t* frags = malloc(sizeof(uint32_t) * n_labels);
  if(n_labels == 0) {
    return frags;
  }
  frags[0] = labels[0].position / bin_size;
  int i;
  for(i = 1; i < n_labels; i++) {
    frags[i] = (labels[i].position - labels[i-1].position) / bin_size;
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
  int i, j;

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  int absent;
  khint_t bin; // hash bin (result of kh_put)
  uint8_t* frags = get_fragments(labels, n_labels, bin_size);
  for(i = 0; i <= n_labels-k; i++) {
    //khint_t qgram = qgram_hash((frags+i), k); // khint_t is probably u32
    khint_t qgram = xratio_hash((labels+i), 100); // khint_t is probably u32
    //printf("# %dth %d-gram: %u\n", i, k, qgram);

    // insert qgram:readId,i into db
    bin = kh_put(qgramHash, db, qgram, &absent);
    if(absent) { // bin is empty (unset)
      kv_init(kh_value(db, bin));
    }
    readPos r;
    r.readNum = (read_id << 1); // forward strand since its padded with a 0
    r.pos = i;
    kv_push(readPos, kh_value(db, bin), r);
  }

  return 0;
}

// if reverse is true (1), forward and reverse signatures will be adjacent such that
// signature of read X forward is at index (2*X), and reverse is (2*X + 1)
void build_hash_db(cmap c, int k, khash_t(qgramHash) *db, int readLimit, int bin_size) {

  uint32_t f = 0;
  while (f < c.n_maps) {

    int res = insert_rmap(c.labels[f], c.map_lengths[f], f, k, 0, db, bin_size); // forward strand only right now
    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
}

khash_t(matchHash)* lookup(label* labels, size_t n_labels, uint32_t read_id, int k, unsigned char reverse, khash_t(qgramHash) *db, int max_qgrams, int bin_size) {
  int i, j, m, absent;

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  khash_t(matchHash) *hits = kh_init(matchHash); // allocate hash table
  khint_t bin; // hash bin (result of kh_put)
  uint8_t* frags = get_fragments(labels, n_labels, bin_size);
  if(n_labels < k) return hits;
  for(i = 0; i <= n_labels-k; i++) {
    //khint_t qgram = qgram_hash((frags+i), k); // khint_t is probably u32
    khint_t qgram = xratio_hash((labels+i), 100); // khint_t is probably u32
    //printf("# %dth %d-gram: %u\n", i, k, qgram);

    bin = kh_get(qgramHash, db, qgram);
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
      posPair ppair; // to store matching query/target positions
      ppair.qpos = i;
      ppair.tpos = kv_A(matches, m).pos;
      kv_push(posPair, kh_value(hits, bin), ppair);
    }
  }

  return hits;
}

void query_db(cmap b, int k, khash_t(qgramHash) *db, cmap c, int readLimit, int max_qgrams, int threshold, int bin_size) {
  int i, j;
  uint32_t target;
  int max_chains = 10; // this can be a parameter
  int match_score = 4; // ? idk what to do with this
  int max_gap = 10; // need to test/refine this
  int min_chain_length = 3; // need to test/refine this

  uint32_t f = 0;
  while (f < b.n_maps) {
    unsigned char qrev = 0; // TODO: forward strand only right now - the database is also currently only fw, eventually we should query both

    printf("# Hashing fragment of size %d with %d nicks\n", b.ref_lengths[f], b.map_lengths[f]);
    khash_t(matchHash) *hits = lookup(b.labels[f], b.map_lengths[f], f, k, qrev, db, max_qgrams, bin_size); // forward strand only right now

    // assess hits for each target
    pairVec target_hits;
    khint_t iter;
    // from the definition of kh_foreach
    for (iter = kh_begin(hits); iter != kh_end(hits); ++iter) {
      if (!kh_exist(hits, iter)) continue;
      target = kh_key(hits, iter);
      target_hits = kh_val(hits, iter);

      //sort anchor pairs by target pos increasing
      ks_mergesort(pos_pair_cmp, kv_size(target_hits), target_hits.a, 0);

      pairVec* chains = chain(&target_hits, max_chains, match_score, max_gap, min_chain_length);
      int n_chains = 0;
      // count chains (if fewer than max_chains, the chains array is terminated by an empty chain vector)
      for(i = 0; i < max_chains; i++) {
        if(kv_size(chains[i]) > 0)
          n_chains = i+1;
        else
          break;
      }
      printf("hit target %d %d times\n", target, kv_size(target_hits));

      printf("%d chains found with anchor sizes: ", n_chains);
      for(i = 0; i < n_chains; i++) {
        printf("%d, ", kv_size(chains[i]));
      }
      printf("\n");

      //if(kv_size(target_hits) >= threshold && target != f) {

      for(j = 0; j < n_chains; j++) {

        // dynamic time warping
        // extract ref labels - expand bounds to encompass unmatched labels within query range
        int rst = kv_A(chains[j], 0).tpos;
        // esimated start position on ref is (anchor[0]_ref_pos - anchor[0]_query_pos)
        int est_rst = c.labels[target][kv_A(chains[j], 0).tpos].position - b.labels[f][kv_A(chains[j], 0).qpos].position;
        while(rst > 0 && c.labels[target][rst].position > est_rst)
          rst--;
        int ren = kv_A(chains[j], kv_size(chains[j])-1).tpos;
        // estimated end position on ref is (anchor[n]_ref_pos + (query_length - anchor[n]_query_pos))
        int est_ren = c.labels[target][kv_A(chains[j], kv_size(chains[j])-1).tpos].position + (b.labels[f][b.map_lengths[f]-1].position - b.labels[f][kv_A(chains[j], kv_size(chains[j])-1).qpos].position);
        while(ren < c.map_lengths[target]-1 && c.labels[target][ren].position < est_ren)
          ren++;
        //fprintf(stderr, "est ref pos %u - %u\n", est_rst, est_ren);
        //fprintf(stderr, "r indices %u - %u\n", rst, ren);

        pathvec path;
        kv_init(path);
        // get fragment distances for DTW (no dicretization)
        uint32_t* qfrags = u32_get_fragments(b.labels[f], b.map_lengths[f], 1);
        uint32_t* rfrags = u32_get_fragments(c.labels[target]+rst, ren-rst+1, 1);
        result aln = dtw(qfrags, rfrags, b.map_lengths[f], ren-rst+1, &path, -1, -1, 1000); // ins_score, del_score, neutral_deviation
        free(qfrags);
        free(rfrags);
        if(aln.failed) {
          printf("q %d : t %d DTW failed\n", f, target);
          continue;
        }

        if(!aln.failed) { // && aln.score >= threshold) {
          // q, t, qstart, qend, qlen, tstart, tend, tlen, score
          printf("DTW result: %d,%d,%d,%d,%d,%d,%d,%d,%f\n", f, target, aln.qstart, aln.qend, b.map_lengths[f], aln.tstart, aln.tend, ren-rst+1, aln.score);
          printf("path: ");
          for(i = 0; i < kv_size(path); i++) {
            printf("%c ", kv_A(path, i) == 0 ? '.' : (kv_A(path, i) == 1 ? 'I' : 'D'));
          }
          printf("\n");
        }

        printf("%d,%d,%d,%d", f, qrev, target, kv_size(chains[j]));
        for(i = 0; i < kv_size(chains[j]); i++)
          printf(",%u(%u):%u(%u)", kv_A(chains[j], i).qpos, b.labels[f][kv_A(chains[j], i).qpos].position, kv_A(chains[j], i).tpos, c.labels[target][kv_A(chains[j], i).tpos].position);
        printf("\n");
      }
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
int hash_cmap(cmap b, cmap c, int q, int threshold, int max_qgrams, int readLimit, int bin_size) {

  // ------------------------- Create hash database -----------------------------

  time_t t0 = time(NULL);

  // read BNX file, construct hash, including only forward direction
  // -- maybe assess doing this the opposite way at a later date, I'm not sure which will be faster
  khash_t(qgramHash) *db = kh_init(qgramHash);

  printf("# Hashing %d cmap fragments\n", c.n_maps);
  build_hash_db(c, q, db, readLimit, bin_size);

  time_t t1 = time(NULL);
  printf("# Hashed rmaps in %d seconds\n", (t1-t0));
  t0 = t1;
  // -------------------------------------------------------------------------------

  // ---------------------------- Look up queries in db ------------------------------
  printf("# Querying %d bnx fragments\n", b.n_maps);
  query_db(b, q, db, c, readLimit, max_qgrams, threshold, bin_size);

  t1 = time(NULL);
  printf("# Queried and output in %d seconds\n", (t1-t0));
  // ----------------------------------------------------------------------------------------

  // TODO: clean up vectors in each hash bin
  kh_destroy(qgramHash, db);
  return 0;
}
