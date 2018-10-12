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
#include "lsh.h"
#include "hash.h"


/*
 * nicks: an array of nick positions as produced by bntools (bn_file.c)
 * k: q-gram size
 * h: number of unique hashes to compute
 * denom: number of divide every raw value by - any results that exceed 255 will be clipped to 255
 *   - for denom = 1000, any fragments over 255Kb will appear to be 255
 * hash_seeds: the random seed set to create unique hashes (|hash_seeds| == h)
 * reverse: 0 - no, 1 - yes
 *
 * returns: an array of Min values (length is h)
 */
Min* minhash (byteVec frags, int k, int h, uint32_t denom, uint32_t hash_seeds[], unsigned char reverse) {

  int i, j;
  int slen = kv_size(frags);

  // initialize minimums, one per hash seed
  Min* m = (Min*)malloc(sizeof(Min)*h);

  for(j = 0; j < h; j++) {
    // initialize [if a whole read is Gs, it will actually have a minhash value == UINT32_MAX, in which case the position will be 0]
    m[j].hash = UINT32_MAX;
    m[j].pos = 0;
  }

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  for(i = 0; i <= slen-k; i++) {
    khint_t qgram = qgram_hash((frags.a+i), k, -1); // khint_t is probably u32
    //printf("# %dth %d-gram: %d\n", i, k, qgram);

    for(j = 0; j < h; j++) {
      khint_t xor_qgram = qgram ^ hash_seeds[j];
      if(xor_qgram < m[j].hash) {
        m[j].hash = xor_qgram;
        m[j].pos = i;
      }
    }
  }

  return m;
}

// if reverse is true (1), forward and reverse signatures will be adjacent such that
// signature of read X forward is at index (2*X), and reverse is (2*X + 1)
void hash_signatures(byteVec *frags, size_t n_frags, int k, int h, uint32_t hash_seeds[], khash_t(qgramHash) **min_db, minVec *min_queries, int readLimit) {

  int l, i, ret_val;
  uint32_t denom = 100; // bin fragments by 100bp

  khint_t bin; // hash bin (result of kh_put)

  uint32_t f = 0;
  while (f < n_frags) {

    Min *m = minhash(frags[f], k, h, denom, hash_seeds, 0);

    // only add hashes to the lookup dictionary if |sub-frags| >= k, otherwise they would be UINT32_MAX
    if(min_db != NULL && kv_size(frags[f]) >= k) {
      for(i = 0; i < h; i++) {
        bin = kh_put(qgramHash, min_db[i], m[i].hash, &ret_val);
        if(ret_val == 1) { // bin is empty (unset)
          kv_init(kh_value(min_db[i], bin)); // kh_value should already be defined as type matchVec
        }
        readPos r;
        r.readNum = (f << 1);
        r.pos = m[i].pos;
        kv_push(readPos, kh_value(min_db[i], bin), r);
      }
    }
    if(min_queries != NULL) {
      kv_push(Min*, *min_queries, m);
      m = minhash(frags[f], k, h, denom, hash_seeds, 1);
      kv_push(Min*, *min_queries, m);
    } else {
      free(m);
    }
    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
}

int free_lsh(khash_t(qgramHash) **min_db, int h, minVec *min_queries) {
  int i;
  for(i = 0; i < h; i++)
    kh_destroy(qgramHash, min_db[i]);
  kv_destroy(*min_queries);
  return 0;
}

int ovl_rmap(byteVec *frags, size_t n_frags, int q, int h, int seed, int threshold, int max_qgrams, int readLimit) {

  // construct array of hash seeds
  srand(seed); // seed random number generator
  // RAND_MAX must be at least 2**15, mine is 2**31
  // as long as we don't duplicate seeds, it should be fine
  uint32_t hash_seeds[h];
  int i;
  for(i = 0; i < h; i++) {
    hash_seeds[i] = rand(); // 0 to RAND_MAX
  }

  // ------------------------- Create hash database -----------------------------

  time_t t0 = time(NULL);

  // read BNX file, construct hash, including only forward direction
  // -- maybe assess doing this the opposite way at a later date, I'm not sure which will be faster
  khash_t(qgramHash) *min_db[h];
  for(i = 0; i < h; i++) {
    min_db[i] = kh_init(qgramHash); // allocate hash table
  }
  minVec min_queries;
  kv_init(min_queries);

  printf("# Hashing %d rmap fragments\n", n_frags);
  hash_signatures(frags, n_frags, q, h, hash_seeds, min_db, &min_queries, readLimit); // load both db and queires from single bnx

  time_t t1 = time(NULL);
  printf("# Hashed %d rmaps in %d seconds\n", min_queries.n/2, (t1-t0));
  t0 = t1;
  // -------------------------------------------------------------------------------

  // ---------------------------- Look up queries in db ------------------------------
  int j, qidx, qrev, m, ret_val;
  uint32_t target;
  khint_t bin; // hash bin (result of kh_put, kh_get)
  Min *qmin;
  matchVec matches;
  khash_t(matchHash) *overlaps;
  for(i = 0; i < min_queries.n/2; i++) {
    for(qrev = 0; qrev <= 1; qrev++) {
      if(qrev == 1) {
        qidx = i*2 + 1;
      } else {
        qidx = i*2;
      }

      // this is the minhash sketch for qidx
      // - if the query was smaller than k, all of the min values will be UINT32_MAX
      //   and it won't hit anything because we didn't put empties in the lookup dict
      qmin = (Min*)(kv_A(min_queries, qidx));

      overlaps = kh_init(matchHash); // allocate hash table
      for(j = 0; j < h; j++) {
        bin = kh_get(qgramHash, min_db[j], qmin[j].hash);
        if(bin == kh_end(min_db[j])) // key not found, IDK what happens if you don't test this
          continue;
        matches = kh_val(min_db[j], bin);
        if(matches.n > max_qgrams) { // repetitive, ignore it
          continue;
        }
        for(m = 0; m < matches.n; m++) {
          bin = kh_put(matchHash, overlaps, matches.a[m].readNum>>1, &ret_val); // >>1 removes the fw/rv bit, which is always fw(0) right now
          if(ret_val == 1) { // bin is empty (unset)
            kv_init(kh_value(overlaps, bin)); // kh_value should already be defined as type matchVec
          }
          posPair ppair; // to store matching query/target positions
          ppair.qpos = qmin[j].pos;
          ppair.tpos = matches.a[m].pos;
          kv_push(posPair, kh_value(overlaps, bin), ppair);
        }
      }
      // count hits and compute offset for each target
      pairVec offsets;
      khint_t iter;
      // from the definition of kh_foreach
      for (iter = kh_begin(overlaps); iter != kh_end(overlaps); ++iter) {
        if (!kh_exist(overlaps, iter)) continue;
        target = kh_key(overlaps, iter);
        offsets = kh_val(overlaps, iter);
        if(offsets.n >= threshold && target != i) {
          /*
           * To compute offset directly:
           *
          offset = 0;
          for(i = 0; i < offsets.n; i++) {
            offset = offset + (offset.a[i].tpos - offsets.a[i].qpos);
          }
          offset = offset / (int)offsets.n; // offsets.n is a size_t, does not play nicely
           *
           */

          printf("%d,%d,%d,%d", i, qrev, target, offsets.n);
          for(j = 0; j < offsets.n; j++)
            printf(",%d:%d", offsets.a[j].qpos, offsets.a[j].tpos);
          printf("\n");
        }
        kv_destroy(offsets); // this frees each value vector after its been evaluated
      }
      kh_destroy(matchHash, overlaps); // this frees the encompassing hash when we're done
    }
  }

  t1 = time(NULL);
  printf("# Compared and output in %d seconds\n", (t1-t0));
  // ----------------------------------------------------------------------------------------

  free_lsh(min_db, h, &min_queries);
  return 0;
}
