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
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

#include "cmap.h"
#include "bnx.h"
//#include "hash.h"
//#include "lsh.h"
//#include "dtw.h"
#include "sim.h"
#include "digest.h"

void usage() {
  printf("Usage: rekit [command] [options]\n");
  printf("Commands:\n");
  printf("  ovl: compute MinHash/pairwise Jaccard similarity\n");
  printf("  aln: compute dynamic time warping glocal (overlap) alignments\n");
  printf("  hsh: compute full q-gram intersection using a hash table\n");
  printf("  sim: simulate molecules\n");
  printf("  dig: in silico digestion\n");
  printf("Options:\n");
  printf("  ovl -b\n");
  printf("  aln -b\n");
  printf("  hsh -b\n");
  printf("  sim -frx --break-rate --fn --fp --min-frag --stretch-mean --stretch-std\n");
  printf("  dig -fr\n");
  printf("    -b: bnx: A single BNX file containing molecules\n");
  printf("    -c: cmap: A single CMAP file\n");
  printf("    -f: fasta: Reference sequence to simulate from\n");
  printf("    -r: cutseq: Recognition/label site sequence\n");
  printf("    -q: Size of q-gram/k-mer to hash\n");
  printf("    -h: Number of hash functions to apply\n");
  printf("    -s: Seed to random number generator\n");
  printf("    -t: threshold: Minimum number of q-grams to declare a match\n");
  printf("    -m: max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored\n");
  printf("    -d: threshold: Score threshold to report alignment\n");
  printf("    -x: Simulated molecule coverage\n");
  printf("  sim options:\n");
  printf("    --break-rate: Probability of genome fragmentation per locus (default: 0.000005)\n");
  printf("    --fn: Probability of missed label at true restriction site (default: 0.1)\n");
  printf("    --fp: Probability of false-positive label (default: 0.05)\n");
  printf("    --stretch-mean: Fragment stretch mean (default: 1.0)\n");
  printf("    --stretch-std: Fragment stretch standard deviation (default: 0.05)\n");
  printf("    --min-frag: Minimum detectable fragment size\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "break-rate",      required_argument, 0, 0 },
  { "fn",              required_argument, 0, 0 },
  { "fp",              required_argument, 0, 0 },
  { "min-frag",        required_argument, 0, 0 },
  { "stretch-mean",    required_argument, 0, 0 },
  { "stretch-std",     required_argument, 0, 0 },
  { 0, 0, 0, 0}
};

int main(int argc, char *argv[]) {

  char* bnx_file = NULL; // .bnx file path/name
  char* fasta_file = NULL; // .fasta file path/name
  char* cmap_file = NULL; // .cmap file path/name (either reference or consensus)
  char* restriction_seq = NULL; // restriction enzyme or label recognition sequence (must also be reverse complemented if not symmetrical)
  int q = 5; // q-gram size
  int h = 10; // number of hashes
  int verbose = 0;
  int threshold = 0;
  int seed = 0; // made this up
  int max_qgrams = 1000; // made this up

  int coverage = 0;
  float break_rate = 0.000005; // one every 200Kb
  float fn = 0.1; // 10% false negative (missed) labels
  float fp = 0.05; // 5% false positive (extra) labels
  float min_frag = 500; // minimum reported gap between labels
  float stretch_mean = 1.0; // no uniform stretch
  float stretch_std = 0.05; // 5% random (normally distributed) sizing error

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "b:c:q:h:f:r:t:m:vx:", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'b':
        bnx_file = optarg;
        break;
      case 'c':
        cmap_file = optarg;
        break;
      case 'q':
        q = atoi(optarg);
        break;
      case 'h':
        h = atoi(optarg);
        break;
      case 'f':
        fasta_file = optarg;
        break;
      case 'r':
        restriction_seq = optarg;
        break;
      case 't':
        threshold = atoi(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'x':
        coverage = atoi(optarg);
        break;
      case '?':
        if (optopt == 'b' || optopt == 'c' || optopt == 'q' || optopt == 'h' || optopt == 'r' || optopt == 'f' || optopt == 't' || optopt == 'm' || optopt == 'x')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        if (long_idx == 0) break_rate = atof(optarg); // --break-rate
        else if (long_idx == 1) fn = atof(optarg); // --fn
        else if (long_idx == 2) fp = atof(optarg); // --fp
        else if (long_idx == 3) min_frag = atoi(optarg); // --min-frag
        else if (long_idx == 4) stretch_mean = atof(optarg); // --stretch-mean
        else if (long_idx == 5) stretch_std = atof(optarg); // --stretch-std
      default:
        usage();
        return 1;
    }
  }

  // CMAP test
  //cmap c = read_cmap(cmap_file);
  //printf("CMAP '%s': %d maps w/%d recognition sites\n", cmap_file, c.n_maps, c.n_rec_seqs);
  //write_cmap(&c, fopen("test.cmap", "w"));

  int index;
  char* command = NULL;
  for (index = optind; index < argc; index++) {
    if(index == optind) {
      command = argv[index];
    }
  }
  if(command == NULL) {
    fprintf(stderr, "Missing command (sim, dig, ovl, aln, hsh)\n");
    return 1;
  }

  cmap c;
  size_t n_frags;

  int i = 0;
  int ret = 1;

  if(strcmp(command, "dig") == 0) {
    fprintf(stderr, "-- Running in silico digest --\n");
    if(fasta_file == NULL) {
      fprintf(stderr, "FASTA file required (-f)\n");
      return 1;
    }
    if(restriction_seq == NULL) {
      fprintf(stderr, "Restriction sequence is required (-r)\n");
      return 1;
    }
    // make a list of restriction seqs - that's what digest wants
    char** rseqs = malloc(1 * sizeof(char*));
    rseqs[0] = restriction_seq;
    c = digest_fasta(fasta_file, rseqs, 1);
    ret = write_cmap(&c, stdout);
  }

  if(strcmp(command, "ovl") == 0) {
    printf("# Loading '%s'...\n", bnx_file);
    time_t t0 = time(NULL);
    c = read_bnx(bnx_file);

    time_t t1 = time(NULL);
    printf("# Loaded %d molecules in %d seconds\n", c.n_maps, (t1-t0));

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    //byteVec *frags = nicks_to_frags_bin(&map, 1000, -1, 1, 1, 1000);

    //int ret = ovl_rmap(frags, n_frags, q, h, seed, threshold, max_qgrams, -1);

    // free
    /*
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
    */
  }

  else if(strcmp(command, "aln") == 0) {
    if(argc < 4) {
      printf("Not enough arguments for 'aln'.\n\n");
      usage();
      return 1;
    }

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    /*
    u32Vec *frags = nicks_to_frags(&map, -1, 1, 1, 1000);

    int ret = dtw_rmap(frags, n_frags, threshold);

    // free
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
    */
  }

  else if(strcmp(command, "hsh") == 0) {
    if(argc < 6) {
      printf("Not enough arguments for 'aln'.\n\n");
      usage();
      return 1;
    }

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    /*
    byteVec *frags = nicks_to_frags_bin(&map, 1000, -1, 1, 1, 1000);

    int ret = hash_rmap(frags, n_frags, q, threshold, max_qgrams, -1);

    // free
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
    */
  }

  else if(strcmp(command, "sim") == 0) {
    if(fasta_file == NULL) {
      fprintf(stderr, "FASTA file required (-f)\n");
      return 1;
    }
    if(restriction_seq == NULL) {
      fprintf(stderr, "Restriction sequence is required (-r)\n");
      return 1;
    }
    if(coverage == 0) {
      fprintf(stderr, "Coverage is required (-x)\n");
      return 1;
    }
    // make a list of restriction seqs - that's what digest wants
    char** rseqs = malloc(1 * sizeof(char*));
    rseqs[0] = restriction_seq;

    fprintf(stderr, "-- Running optical mapping simulation --n");
    c = simulate_bnx(fasta_file, rseqs, 1, break_rate, fn, fp, stretch_mean, stretch_std, min_frag, coverage);
    ret = write_cmap(&c, stdout);
  }

  // free everything
  return ret;
}
