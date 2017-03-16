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
#include <time.h>

#include "bnx.h"
#include "rmap.h"
#include "hash.h"
#include "lsh.h"
#include "dtw.h"

void usage() {
  printf("Usage: rekit [command] [options]\n");
  printf("Commands:\n");
  printf("  ovl: compute MinHash/pairwise Jaccard similarity\n");
  printf("  aln: compute dynamic time warping glocal (overlap) alignments\n");
  printf("  hsh: compute full q-gram intersection using a hash table\n");
  printf("Options:\n");
  printf("  ovl <bnx> <q> <h> <seed> <threshold> <max_qgram_hits>\n");
  printf("    bnx: A single BNX file containing rmaps\n");
  printf("    q: Size of q-gram/k-mer to hash\n");
  printf("    h: Number of hash functions to apply\n");
  printf("    seed: Seed to random number generator\n");
  printf("    threshold: Minimum number of q-grams to declare a match\n");
  printf("    max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored\n");
  printf("  aln <bnx> <threshold>\n");
  printf("    bnx: A single BNX file containing rmaps\n");
  printf("    threshold: Score threshold to report alignment\n");
  printf("  hsh <bnx> <q> <threshold> <max_qgram_hits>\n");
  printf("    bnx: A single BNX file containing rmaps\n");
  printf("    q: Size of q-gram/k-mer to hash\n");
  printf("    threshold: Minimum number of q-grams to declare a match\n");
  printf("    max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored\n");
}

/*
 * Takes a set of rmaps (from one BNX file)
 * and returns an array of byteVecs, one for each fragment
 * where the values are discretized, potentially trimmed, sub-fragment sizes (between nick sites)
 */
byteVec* nicks_to_frags_bin(rmap *map, uint32_t denom, int readLimit, int ltrim, int rtrim) {
  byteVec *frags = (byteVec*)malloc(sizeof(byteVec) * kv_size(map->fragments));
  int i;
  uint32_t f = 0;
  while (f < kv_size(map->fragments)) {

    nickVec nicks = kv_A(map->fragments, f).nicks;
    int slen = kv_size(nicks);

    // build ordered kvec of fragment sizes instead of nick positions
    // and bin their values into nbins
    kv_init(frags[f]);
    for(i = ltrim; i < slen-rtrim; i++) {
      uint32_t frg_size = (kv_A(nicks, i).pos - (i>0 ? kv_A(nicks, i-1).pos : 0));
      uint8_t bin = (frg_size / denom > 255) ? 255 : (frg_size / denom);
      kv_push(uint8_t, frags[f], bin);
    }

    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
  return frags;
}

u32Vec* nicks_to_frags(rmap *map, int readLimit, int ltrim, int rtrim) {
  u32Vec *frags = (u32Vec*)malloc(sizeof(u32Vec) * kv_size(map->fragments));
  int i;
  uint32_t f = 0;
  while (f < kv_size(map->fragments)) {

    nickVec nicks = kv_A(map->fragments, f).nicks;
    int slen = kv_size(nicks);

    // build ordered kvec of fragment sizes instead of nick positions
    // and bin their values into nbins
    kv_init(frags[f]);
    for(i = ltrim; i < slen-rtrim; i++) {
      uint32_t frg_size = (kv_A(nicks, i).pos - (i>0 ? kv_A(nicks, i-1).pos : 0));
      kv_push(uint32_t, frags[f], frg_size);
    }

    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
  return frags;
}


int main(int argc, char *argv[]) {
  if(argc < 3) {
    printf("Not enough arguments - specify a command and BNX file first.\n\n");
    usage();
    return -1;
  }
  char *command = argv[1];
  char *bnx_file = argv[2];

  // test if files exist
  FILE* fp;
  fp = fopen(bnx_file, "r");
  if(fp == NULL) {
    printf("File '%s' does not exist\n", bnx_file);
    return -1;
  }
  fclose(fp);

  printf("# Loading '%s' into rmap\n", bnx_file);
  time_t t0 = time(NULL);
  rmap map = bn_load(bnx_file);

  size_t n_frags = kv_size(map.fragments);
  time_t t1 = time(NULL);
  printf("# Loaded in %d seconds\n", (t1-t0));

  int i = 0;
  int ret = 1;

  if(strcmp(command, "ovl") == 0) {
    if(argc < 8) {
      printf("Not enough arguments for 'ovl'.\n\n");
      usage();
      return -1;
    }
    int q = atoi(argv[3]);
    int h = atoi(argv[4]);
    int seed = atoi(argv[5]);
    int threshold = atoi(argv[6]);
    int max_qgrams = atoi(argv[7]);

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    byteVec *frags = nicks_to_frags_bin(&map, 1000, -1, 1, 1);

    int ret = ovl_rmap(frags, n_frags, q, h, seed, threshold, max_qgrams, -1);

    // free
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
  }

  else if(strcmp(command, "aln") == 0) {
    if(argc < 4) {
      printf("Not enough arguments for 'aln'.\n\n");
      usage();
      return -1;
    }
    int threshold = atoi(argv[3]);

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    u32Vec *frags = nicks_to_frags(&map, -1, 1, 1);

    int ret = dtw_rmap(frags, n_frags, threshold);

    // free
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
  }

  else if(strcmp(command, "hsh") == 0) {
    if(argc < 6) {
      printf("Not enough arguments for 'aln'.\n\n");
      usage();
      return -1;
    }
    int q = atoi(argv[3]);
    int threshold = atoi(argv[4]);
    int max_qgrams = atoi(argv[5]);

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    byteVec *frags = nicks_to_frags_bin(&map, 1000, -1, 1, 1);

    int ret = hash_rmap(frags, n_frags, q, threshold, max_qgrams, -1);

    // free
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
  }

  rmap_free(&map);
  return ret;
}
