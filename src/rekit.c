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
  time_t t1 = time(NULL);
  printf("# Loaded in %d seconds\n", (t1-t0));

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

    int ret = ovl_rmap(map, q, h, seed, threshold, max_qgrams, -1);
  }
  else if(strcmp(command, "aln") == 0) {
    if(argc < 4) {
      printf("Not enough arguments for 'aln'.\n\n");
      usage();
      return -1;
    }
    int threshold = atoi(argv[3]);

    int ret = dtw_rmap(map, threshold);
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

    int ret = hash_rmap(map, q, threshold, max_qgrams, -1);
  }

  return ret;
}
