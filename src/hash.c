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
#include "rmap.h"
#include "hash.h"


/*
 * nicks: an array of nick positions as produced by bntools (bn_file.c)
 * k: q-gram size
 * reverse: 0 - no, 1 - yes
 *
 * returns: 0 if successful, else 1
 */
int insert_rmap(byteVec frags, uint32_t read_id, int k, unsigned char reverse, khash_t(qgramHash) *db) {
  int i, j;
  int slen = kv_size(frags);

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  int ret_val;
  khint_t bin; // hash bin (result of kh_put)
  for(i = 0; i <= slen-k; i++) {
    khint_t qgram = qgram_hash((frags.a+i), k); // khint_t is probably u32
    //printf("# %dth %d-gram: %d\n", i, k, qgram);

    // insert qgram:readId,i into db
    bin = kh_put(qgramHash, db, qgram, &ret_val);
    if(ret_val == 1) { // bin is empty (unset)
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
void build_hash_db(byteVec *frags, size_t n_frags, int k, khash_t(qgramHash) *db, int readLimit) {

  uint32_t f = 0;
  while (f < n_frags) {

    int res = insert_rmap(frags[f], f, k, 0, db); // forward strand only right now
    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
}

khash_t(matchHash)* lookup(byteVec frags, uint32_t read_id, int k, unsigned char reverse, khash_t(qgramHash) *db, int max_qgrams) {
  int i, j, m, ret_val;
  int slen = kv_size(frags);

  // nick pos values are ints, rounded from the double in the bnx file, and may be 4 or 8 bytes long
  // they are related to the length of the whole fragment, so typically max out in the 100s of thousands (avg ~200k)
  // there is always a nick value given for the END of the fragment, but not one at 0

  khash_t(matchHash) *hits = kh_init(matchHash); // allocate hash table
  khint_t bin; // hash bin (result of kh_put)
  for(i = 0; i <= slen-k; i++) {
    khint_t qgram = qgram_hash((frags.a+i), k); // khint_t is probably u32
    //printf("# %dth %d-gram: %d\n", i, k, qgram);

    bin = kh_get(qgramHash, db, qgram);
    if(bin == kh_end(db)) // key not found, IDK what happens if you don't test this
      continue;
    matchVec matches = kh_val(db, bin);
    if(kv_size(matches) > max_qgrams) { // repetitive, ignore it
      continue;
    }
    for(m = 0; m < kv_size(matches); m++) {
      bin = kh_put(matchHash, hits, kv_A(matches, m).readNum>>1, &ret_val); // >>1 removes the fw/rv bit, which is always fw(0) right now
      if(ret_val == 1) { // bin is empty (unset)
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

void query_db(byteVec *frags, size_t n_frags, int k, khash_t(qgramHash) *db, int readLimit, int max_qgrams, int threshold) {
  int i;
  uint32_t target;

  uint32_t f = 0;
  while (f < n_frags) {
    unsigned char qrev = 0; // forward strand only right now - the database is also currently only fw, eventually we should query both

    //printf("# Hashing fragment of size %d with %d nicks\n", kv_A(map.fragments, f).size, kv_size(kv_A(map.fragments, f).nicks));
    khash_t(matchHash) *hits = lookup(frags[f], f, k, qrev, db, max_qgrams); // forward strand only right now

    // assess hits for each target
    pairVec target_hits;
    khint_t iter;
    // from the definition of kh_foreach
    for (iter = kh_begin(hits); iter != kh_end(hits); ++iter) {
      if (!kh_exist(hits, iter)) continue;
      target = kh_key(hits, iter);
      target_hits = kh_val(hits, iter);
      if(kv_size(target_hits) >= threshold && target != f) {

        printf("%d,%d,%d,%d", f, qrev, target, kv_size(target_hits));
        for(i = 0; i < kv_size(target_hits); i++)
          printf(",%d:%d", kv_A(target_hits, i).qpos, kv_A(target_hits, i).tpos);
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
int hash_rmap(byteVec *frags, size_t n_frags, int q, int threshold, int max_qgrams, int readLimit) {

  // ------------------------- Create hash database -----------------------------

  time_t t0 = time(NULL);

  // read BNX file, construct hash, including only forward direction
  // -- maybe assess doing this the opposite way at a later date, I'm not sure which will be faster
  khash_t(qgramHash) *db = kh_init(qgramHash);

  printf("# Hashing %d rmap fragments\n", n_frags);
  build_hash_db(frags, n_frags, q, db, readLimit);

  time_t t1 = time(NULL);
  printf("# Hashed rmaps in %d seconds\n", (t1-t0));
  t0 = t1;
  // -------------------------------------------------------------------------------

  // ---------------------------- Look up queries in db ------------------------------
  printf("# Querying %d rmap fragments\n", n_frags);
  query_db(frags, n_frags, q, db, readLimit, max_qgrams, threshold);

  t1 = time(NULL);
  printf("# Queried and output in %d seconds\n", (t1-t0));
  // ----------------------------------------------------------------------------------------

  // TODO: clean up vectors in each hash bin
  kh_destroy(qgramHash, db);
  return 0;
}
