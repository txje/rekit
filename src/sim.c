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

/*
 * A tool to simulate Bionano restriction mapping of an input genome (fasta)
 * as realistically as possible
 *
 * This seems like a good place to start: https://github.com/pingchen09990102/BMSIM
 *
 * Their procedure goes as follows:
 * 1. Probabilistic fragmentation at marked "fragile" sites
 * 2. Random fragmentation everywhere else as a homogeneous Poisson process
 * 3. in silico digestion of each fragment
 *    - true labels are a Bernoulli event with probability of success p
 *    - false labels are added as a Poisson process
 * 4. Use stretch scale factor to report sizes of each digested subfragment
 * 5. Chimeras (bi, tri, or quad) are constructed by randomly appending fragments with some probability
 * 6. Use some Gaussian resolution model to decide if close (<~1kb) cut sites are differetiable
 * 7. Assign signal scores to labels? Randomly? (I'm not sure this is strictly necessary)
 * 8. Randomly sample from a huge collection of fragments to achieve the desired coverage
 */

#include <math.h>
#include <float.h>
#include <time.h>
#include "klib/kseq.h"
#include "klib/kvec.h"
#include "klib/kstring.h"
#include "rmap.h"
#include "digest.h"
#include "sim.h"


// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// init kseq struct
KSEQ_INIT(FILE*, fileread);


void fragment_seq(kstring_t* seq, seqVec* frag_seqs, float frag_prob) {
  int st = 0, i;
  for(i = 1; i < seq->l; i++) {
    if(rand() * (1.0 / RAND_MAX) < frag_prob) {
      // we'll actually just return each fragment string as a pointer and length into the reference sequence - remember this.
      kstring_t* frag_seq = malloc(sizeof(kstring_t));
      frag_seq->l = frag_seq->m = i-st;
      frag_seq->s = seq->s+st;
      kv_push(kstring_t*, *frag_seqs, frag_seq);
      st = i;
    }
  }
  // last section
  kstring_t* frag_seq = malloc(sizeof(kstring_t));
  frag_seq->l = frag_seq->m = i-st;
  frag_seq->s = seq->s+st;
  kv_push(kstring_t*, *frag_seqs, frag_seq);
}


// from https://en.wikipedia.org/wiki/Box-Muller_transform
float normal(float sigma, float mu) {
  const float two_pi = 2.0 * 3.14159265358979323;
  static float z0, z1;

  static uint8_t generate = 0;
  generate = !generate;

  if (!generate)
     return z1 * sigma + mu;

  float u1, u2;
  do {
     u1 = rand() * (1.0 / RAND_MAX);
     u2 = rand() * (1.0 / RAND_MAX);
  } while ( u1 <= FLT_EPSILON );

  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}


void bn_map(seqVec seqs, fragVec* frags, char **motifs, size_t n_motifs, float nick_prob, float shear_prob, float stretch_mean, float stretch_std, uint32_t resolution_min) {
  int i, j;
  int nlimit = 100; // break strings of Ns at least 100bp long
  for(i = 0; i < kv_size(seqs); i++) {
    kstring_t* seq = kv_A(seqs, i);

    // do digestion + shearing with some probability
    u32Vec frg;
    kv_init(frg);
    int res = digest(seq->s, seq->l, motifs, n_motifs, nick_prob, shear_prob, nlimit, &frg);

    // now perform modifications for stretching and limited resolution
    u32Vec* modfrg = (u32Vec*)malloc(sizeof(u32Vec));
    kv_init(*modfrg);
    for(j = 0; j < kv_size(frg); j++) {
      // apply normally distributed size variation
      uint32_t f = kv_A(frg, j) + normal(stretch_mean, stretch_std);

      // then include only fragments that exceed some minimum size (typically, ~1kb for Bionano)
      if(f >= resolution_min) {
        kv_push(uint32_t, *modfrg, f);
      }
    }

    kv_push(u32Vec*, *frags, modfrg);
    kv_destroy(frg);
  }
}


void simulate_bnx(char* ref_fasta, float frag_prob, float nick_prob, float shear_prob, float stretch_mean, float stretch_std, uint32_t resolution_min, float coverage) {
  char **motifs = (char**)malloc(sizeof(char*));
  motifs[0] = "ACGTGCA";
  size_t n_motifs = 1;

  float bimera_prob = 0.01;
  float trimera_prob = 0.001;
  float quadramera_prob = 0.0001;

  srand(time(NULL));

  FILE* fp;
  kseq_t* seq;
  uint32_t i;
  int j;

  size_t genome_size = 0;

  // vector to hold optical fragments
  fragVec fragments;
  kv_init(fragments);
  seqVec frag_seqs;
  kv_init(frag_seqs);

  fp = fopen(ref_fasta, "r");
  seq = kseq_init(fp);
  printf("Fasta file: %s\n", ref_fasta);
  printf("Fragmenting and digesting up to %fx coverage\n", coverage*10);

  int l;
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    printf("Reading sequence '%s' (%i bp).\n", seq->name.s, l);
    fflush(stdout);
    genome_size = genome_size + seq->seq.l;

    // digest 10x as many genomes as the coverage we want, then we'll sample down
    int c;
    for(c = 0; c < coverage*10; c++) {
      frag_seqs.n = 0; // reset the intermediate fragment sequences each round
      fragment_seq(&seq->seq, &frag_seqs, frag_prob);
      bn_map(frag_seqs, &fragments, motifs, n_motifs, nick_prob, shear_prob, stretch_mean, stretch_std, resolution_min);
    }
  }
  kseq_destroy(seq);
  fclose(fp);
  printf("Produced %d read fragments.\n", kv_size(fragments));
  fflush(stdout);
  
  // sample randomly down to the requested coverage, include chimeras
  // bimera_prob, trimera_prob, quadramera_prob
  printf("Sampling and chimerizing down to %fx coverage\n", coverage);
  fflush(stdout);
  fragVec observed;
  kv_init(observed);
  size_t f0, f1;
  for(i = 0; i < coverage*genome_size; ) {
    float chimera_prob = rand() * (1.0 / RAND_MAX);
    int chimera = (chimera_prob < quadramera_prob ? 4 : (chimera_prob < trimera_prob ? 3 : (chimera_prob < bimera_prob ? 2 : 1)));

    // sample without replacement
    do {
      f0 = (size_t)(rand() * (1.0 / RAND_MAX) * kv_size(fragments));
      printf("f: %d\n", f0);
    } while(kv_A(fragments, f0) == NULL); // this will never slow down too much since we created 10x as many as we'll ever sample
    printf("picked read fragment %d\n", f0);
    fflush(stdout);

    // append as many other fragments as the chimerism calls for
    while(chimera > 1) {
      printf("chimera: %d\n", chimera);
      fflush(stdout);
      do {
        f1 = rand() * (1.0 / RAND_MAX) * kv_size(fragments);
      } while(kv_A(fragments, f1) == NULL); // this will never slow down too much since we created 10x as many as we'll ever sample
      printf("added read fragment %d to chimera\n", f1);
      fflush(stdout);

      kv_extend(uint32_t, *kv_A(fragments, f0), *kv_A(fragments, f1));
      kv_destroy(*kv_A(fragments, f1));
      kv_A(fragments, f1) = NULL;
      chimera--;
    }

    kv_push(u32Vec*, observed, kv_A(fragments, f0));

    // add to our running coverage
    for(j = 0; j < kv_size(*kv_A(fragments, f0)); j++) {
      i = i + kv_A(*kv_A(fragments, f0), j);
    }
    printf("total coverage: %d\n", i);
    fflush(stdout);

    kv_A(fragments, f0) = NULL; // but we can clear the pointer so that it's not sampled again
  }
}
