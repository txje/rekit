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
#include <getopt.h>

#include "bnx.h"
#include "rmap.h"
#include "cmap.h"
#include "hash.h"
#include "lsh.h"
#include "dtw.h"
#include "sim.h"
#include "digest.h"

void usage() {
  printf("Usage: rekit [command] [options]\n");
  printf("Commands:\n");
  printf("  ovl: compute MinHash/pairwise Jaccard similarity\n");
  printf("  aln: compute dynamic time warping glocal (overlap) alignments\n");
  printf("  hsh: compute full q-gram intersection using a hash table\n");
  printf("  sim: simulate rmaps\n");
  printf("  dig: in silico digestion\n");
  printf("Options:\n");
  printf("  ovl -bqht seed max_qgram_hits\n");
  printf("  aln -t score_threshold\n");
  printf("  hsh -bqt <max_qgram_hits>\n");
  printf("  sim -frm <pfrag> <pnick> <pshear> <stretch_mean> <stretch_std> <resolution> <coverage>\n");
  printf("  dig -frm\n");
  printf("    -b: bnx: A single BNX file containing rmaps\n");
  printf("    -c: cmap: A single CMAP file\n");
  printf("    -f: fasta: Reference sequence to simulate from\n");
  printf("    -r: cutseq: Recognition/label site sequence\n");
  printf("    -m: mod: MOD file mapping reference to background\n");
  printf("    -q: Size of q-gram/k-mer to hash\n");
  printf("    -h: Number of hash functions to apply\n");
  printf("    seed: Seed to random number generator\n");
  printf("    -t: threshold: Minimum number of q-grams to declare a match\n");
  printf("    max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored\n");
  printf("    threshold: Score threshold to report alignment\n");
  printf("  sim options:\n");
  printf("    pfrag: Probability of genome fragmentation per locus\n");
  printf("    pnick: Probability of nick at true restriction site\n");
  printf("    pshear: Probability of false-positive nicking (labeling/shearing)\n");
  printf("    stretch_mean: Fragment stretch mean\n");
  printf("    stretch_std: Fragment stretch standard deviation\n");
  printf("    resolution: Minimum detectable fragment size\n");
  printf("    coverage: Target genome coverage\n");
}

/*
 * Takes a set of rmaps (from one BNX file)
 * and returns an array of byteVecs, one for each fragment
 * where the values are discretized, potentially trimmed, sub-fragment sizes (between nick sites)
 */
byteVec* nicks_to_frags_bin(rmap *map, uint32_t denom, int readLimit, int ltrim, int rtrim, uint32_t min) {
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
      if(bin >= min/denom) {
        kv_push(uint8_t, frags[f], bin);
      }
    }

    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
  return frags;
}

u32Vec* nicks_to_frags(rmap *map, int readLimit, int ltrim, int rtrim, uint32_t min) {
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
      if(frg_size >= min) {
        kv_push(uint32_t, frags[f], frg_size);
      }
    }

    f++;

    if(readLimit > 0 && f >= readLimit) {
      break;
    }
  }
  return frags;
}


int main(int argc, char *argv[]) {

  char* bnx_file = NULL; // .bnx file path/name
  char* fasta_file = NULL; // .fasta file path/name
  char* cmap_file = NULL; // .cmap file path/name (either reference or consensus)
  char* restriction_seq = NULL; // restriction enzyme or label recognition sequence (must also be reverse complemented if not symmetrical)
  char* mod = NULL; // MOD file path/name
  int q = 5; // q-gram size
  int h = 10; // number of hashes
  int verbose = 0;
  int threshold = 0;
  int seed = 0; // made this up
  int max_qgrams = 1000; // made this up

  int opt;
  opterr = 0;
  while ((opt = getopt(argc, argv, "b:c:q:h:f:r:t:m:v")) != -1) {
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
      case 'm':
        mod = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case '?':
        if (optopt == 'b' || optopt == 'c' || optopt == 'q' || optopt == 'h' || optopt == 'r' || optopt == 'f' || optopt == 't' || optopt == 'm')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
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
  char* command;
  for (index = optind; index < argc; index++) {
    if(index == optind) {
      command = argv[index];
    }
  }
  printf("command: %s\n", command);

  rmap map;
  size_t n_frags;

  if(strcmp(command, "dig") == 0) {
    printf("-- in silico digest --\n");
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
    cmap c = digest_fasta(fasta_file, rseqs, 1);
    write_cmap(&c, stdout);
  }

  if(strcmp(command, "sim") == 0) {
    // test if files exist
    FILE* fp;
    fp = fopen(bnx_file, "r");

    if(fp == NULL) {
      fprintf(stderr, "File '%s' does not exist\n", bnx_file);
      return 1;
    }
    fclose(fp);

    printf("# Loading '%s' into rmap\n", bnx_file);
    time_t t0 = time(NULL);
    map = bn_load(bnx_file);

    n_frags = kv_size(map.fragments);
    time_t t1 = time(NULL);
    printf("# Loaded in %d seconds\n", (t1-t0));
  }

  int i = 0;
  int ret = 1;

  if(strcmp(command, "ovl") == 0) {
    if(argc < 8) {
      printf("Not enough arguments for 'ovl'.\n\n");
      usage();
      return -1;
    }

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    byteVec *frags = nicks_to_frags_bin(&map, 1000, -1, 1, 1, 1000);

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

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    u32Vec *frags = nicks_to_frags(&map, -1, 1, 1, 1000);

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

    // convert map of nicks to a simple array of byte vectors representing subfragment lengths
    // discretize by 1kb if doing any hashing, do no dicretization if doing DTW
    // -1 can be changed to a positive to number to limit the total number of fragments that are considered
    byteVec *frags = nicks_to_frags_bin(&map, 1000, -1, 1, 1, 1000);

    int ret = hash_rmap(frags, n_frags, q, threshold, max_qgrams, -1);

    // free
    for(i = 0; i < n_frags; i++) {
      kv_destroy(frags[i]);
    }
    free(frags);
  }

  else if(strcmp(command, "sim") == 0) {
    if(argc < 11) {
      printf("Not enough arguments for 'sim'.\n\n");
      usage();
      return -1;
    }
    char *ref_fasta = argv[3];
    float frag_prob = atof(argv[4]);
    float nick_prob = atof(argv[5]);
    float shear_prob = atof(argv[6]);
    float stretch_mean = atof(argv[7]);
    float stretch_std = atof(argv[8]);
    uint32_t resolution = atoi(argv[9]);
    float coverage = atof(argv[10]);

    simulate_bnx(ref_fasta, frag_prob, nick_prob, shear_prob, stretch_mean, stretch_std, resolution, coverage);
  }

  if(strcmp(command, "sim") != 0) {
    rmap_free(&map);
  }
  return ret;
}
