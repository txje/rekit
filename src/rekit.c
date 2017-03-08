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
#include "lsh.h"
#include "dtw.h"


int main(int argc, char *argv[]) {
  if(argc < 7) {
    printf("Usage: rekit <bnx> <q> <h> <seed> <threshold> <max_qgram_hits>\n");
    printf("  bnx: A single BNX file containing rmaps\n");
    printf("  q: Size of q-gram/k-mer to hash\n");
    printf("  h: Number of hash functions to apply\n");
    printf("  seed: Seed to random number generator\n");
    printf("  threshold: Minimum number of k-mers to declare a match\n");
    printf("  max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bnx_file = argv[1];
  int q = atoi(argv[2]);
  int h = atoi(argv[3]);
  int seed = atoi(argv[4]);
  int threshold = atoi(argv[5]);
  int max_qgrams = atoi(argv[6]);

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

  int ret = ovl_rmap(map, q, h, seed, threshold, max_qgrams, -1);
  return ret;
}
