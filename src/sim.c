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
#include <zlib.h>
#include "klib/kvec.h"
#include "klib/ksort.h"
#include "klib/kstring.h"
#include "digest.h"
#include "sim.h"

#ifndef _kseq_
#define _kseq_

#include "klib/kseq.h"

// init kseq struct
KSEQ_INIT(gzFile, gzread)

#endif

KSORT_INIT_GENERIC(uint32_t)


void fragment_seq(kstring_t* seq, seqVec* frag_seqs, posVec* frag_positions, uint32_t ref, float frag_prob) {
  uint32_t st = 0, i;
  for(i = 1; i < seq->l; i++) {
    if(rand() * (1.0 / RAND_MAX) < frag_prob) {
      // we'll actually just return each fragment string as a pointer and length into the reference sequence - remember this.
      kstring_t* frag_seq = malloc(sizeof(kstring_t));
      frag_seq->l = frag_seq->m = i-st;
      frag_seq->s = seq->s+st;
      //printf("added sequence of length %d\n", frag_seq->l);
      kv_push(kstring_t*, *frag_seqs, frag_seq);
      ref_pos rp;
      rp.ref_id = ref;
      rp.pos = st;
      kv_push(ref_pos, *frag_positions, rp);
      st = i;
    }
  }
  // last section
  kstring_t* frag_seq = malloc(sizeof(kstring_t));
  frag_seq->l = frag_seq->m = i-st;
  frag_seq->s = seq->s+st;
  kv_push(kstring_t*, *frag_seqs, frag_seq);
  ref_pos rp;
  rp.ref_id = ref;
  rp.pos = st;
  kv_push(ref_pos, *frag_positions, rp);
}


// from https://en.wikipedia.org/wiki/Box-Muller_transform
float normal(float mu, float sigma) {
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


void bn_map(seqVec seqs, fragVec* frags, char **motifs, size_t n_motifs, float fn_rate, float fp_rate, float stretch_mean, float stretch_std, uint32_t resolution_min) {
  int i, j, k;
  int nlimit = 100; // break strings of Ns at least 100bp long
  for(i = 0; i < kv_size(seqs); i++) {
    kstring_t* seq = kv_A(seqs, i);

    // do digestion + shearing with some probability
    u32Vec positions;
    kv_init(positions);
    int res = digest(seq->s, seq->l, motifs, n_motifs, 1, 0, nlimit, &positions);
    //printf("molecule %d of %d (%d bp) has %d labels\n", i, kv_size(seqs), seq->l, kv_size(positions));

    int fp = round(kv_size(positions) * normal(fp_rate, 0.01));
    u32Vec fp_pos;
    kv_init(fp_pos);
    for(j = 0; j < fp; j++) {
      kv_push(uint32_t, fp_pos, (uint32_t)round((double)rand() / (double)RAND_MAX * kv_A(positions, kv_size(positions)-1)));
    }
    ks_mergesort(uint32_t, kv_size(fp_pos), fp_pos.a, 0); // sort
    k = 0; // index into fp_pos

    // now perform modifications for FN, FP, stretching, and limited resolution
    u32Vec* modpos = (u32Vec*)malloc(sizeof(u32Vec));
    kv_init(*modpos);
    uint32_t last = 0;
    uint32_t last_stretched = 0;
    for(j = 0; j < kv_size(positions); j++) {
      uint32_t val;
      if(k < kv_size(fp_pos) && kv_A(fp_pos, k) < kv_A(positions, j)) {
        val = kv_A(fp_pos, k);
        k++;
        j--;
      } else {
        val = kv_A(positions, j);
      }
      // apply normally distributed size variation
      uint32_t f = last_stretched + (val - last) * normal(stretch_mean, stretch_std);
      last = val;
      last_stretched = f;

      // then include only fragments that exceed some minimum size (typically, ~1kb for Bionano)
      // and fall above FN rate
      if(f >= resolution_min && ((double)rand() / (double)RAND_MAX) > fn_rate) {
        kv_push(uint32_t, *modpos, f);
      }
    }

    kv_push(u32Vec*, *frags, modpos);
    kv_destroy(positions);

    /*
    for(j = 0; j < kv_size(*modpos); j++)
      printf("%d\t", kv_A(*modpos, j));
    printf("\n");
    */
  }
}


cmap simulate_bnx(char* ref_fasta, char** motifs, size_t n_motifs, float frag_prob, float fn, float fp, float stretch_mean, float stretch_std, uint32_t resolution_min, float coverage) {

  float bimera_prob = 0.01;
  float trimera_prob = 0.0001;
  float quadramera_prob = 0.000001;

  srand(time(NULL));

  gzFile f;
  kseq_t* seq;
  uint64_t i;
  int j;

  uint64_t genome_size = 0;

  // vector to hold optical fragments
  fragVec fragments;
  kv_init(fragments);
  posVec frag_positions;
  kv_init(frag_positions);
  seqVec frag_seqs;
  kv_init(frag_seqs);

  f = gzopen(ref_fasta, "r");
  seq = kseq_init(f);
  fprintf(stderr, "Fasta file: %s\n", ref_fasta);
  fprintf(stderr, "Fragmenting and digesting up to %fx coverage\n", coverage*10);

  int l;
  uint32_t ref = 0; // reference seq ID
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //fprintf(stderr, "Reading sequence '%s' (%i bp).\n", seq->name.s, l);
    genome_size = genome_size + seq->seq.l;

    // digest 10x as many genomes as the coverage we want, then we'll sample down
    for(j = 0; j < coverage*10; j++) {
      frag_seqs.n = 0; // reset the intermediate fragment sequences each round
      fragment_seq(&seq->seq, &frag_seqs, &frag_positions, ref, frag_prob);
      bn_map(frag_seqs, &fragments, motifs, n_motifs, fn, fp, stretch_mean, stretch_std, resolution_min);
    }
    ref++;
  }
  kseq_destroy(seq);
  kv_destroy(frag_seqs);
  gzclose(f);

  for(i = 0; i < kv_size(frag_positions); i++)
    fprintf(stderr, "fragment %d position: %u : %u\n", i, kv_A(frag_positions, i).ref_id, kv_A(frag_positions, i).pos);
  
  // sample randomly down to the requested coverage, include chimeras
  // bimera_prob, trimera_prob, quadramera_prob
  fragVec observed;
  kv_init(observed);
  posVec observed_pos;
  kv_init(observed_pos);
  size_t f0, f1;
  uint64_t target_coverage = (uint64_t)((double)coverage * genome_size);
  //fprintf(stderr, "Sampling and chimerizing down to %fx coverage (~%u bp)\n", coverage, target_coverage);
  for(i = 0; i < target_coverage; ) {
    double chimera_prob = (double)rand() / (double)RAND_MAX;
    int chimera = (chimera_prob < quadramera_prob ? 4 : (chimera_prob < trimera_prob ? 3 : (chimera_prob < bimera_prob ? 2 : 1)));

    // sample without replacement
    // this will never slow down too much since we created 10x as many as we'll ever sample
    do {
      f0 = (size_t)(rand() * (1.0 / RAND_MAX) * kv_size(fragments));
    } while(kv_A(fragments, f0) == NULL || kv_size(*kv_A(fragments, f0)) <= 0); // keep trying if this fragment was already sampled or has no labels

    // append as many other fragments as the chimerism calls for
    while(chimera > 1) {
      do {
        f1 = rand() * (1.0 / RAND_MAX) * kv_size(fragments);
      } while(f1 == f0 || kv_A(fragments, f1) == NULL); // this will never slow down too much since we created 10x as many as we'll ever sample

      // shift all of the f1 positions over to append to f0
      uint32_t last = kv_A(*kv_A(fragments, f0), kv_size(*kv_A(fragments, f0))-1);
      for(j = 0; j < kv_size(*kv_A(fragments, f1)); j++) {
        kv_push(uint32_t, *kv_A(fragments, f0), kv_A(*kv_A(fragments, f1), j) + last);
      }

      kv_destroy(*kv_A(fragments, f1));
      kv_A(fragments, f1) = NULL;
      chimera--;
    }

    /*
    printf("%d\t", chimera);
    for(j = 0; j < kv_size(*kv_A(fragments, f0)); j++)
      printf("%d\t", kv_A(*kv_A(fragments, f0), j));
    printf("\n");
    */

    kv_push(u32Vec*, observed, kv_A(fragments, f0));
    kv_push(ref_pos, observed_pos, kv_A(frag_positions, f0)); // does not record the other parts (if any) of a chimera
    fprintf(stderr, "observing fragment %d position: %u : %u\n", f0, kv_A(frag_positions, f0).ref_id, kv_A(frag_positions, f0).pos);

    // add to our running coverage
    size_t f0_size = kv_size(*kv_A(fragments, f0));
    i += kv_A(*kv_A(fragments, f0), kv_size(*kv_A(fragments, f0))-1);

    kv_A(fragments, f0) = NULL; // but we can clear the pointer so that it's not sampled again
  }

  // clean up unsampled fragments
  for(i = 0; i < kv_size(fragments); i++)
    if(kv_A(fragments, i) != NULL)
      kv_destroy(*kv_A(fragments, i));
  kv_destroy(frag_positions);
  kv_destroy(fragments);

  cmap c;
  init_cmap(&c);
  c.n_rec_seqs = n_motifs;
  c.rec_seqs = motifs;
  for(i = 0; i < kv_size(observed); i++) {
    add_map(&c, kv_A(observed, i)->a, kv_size(*kv_A(observed, i)), 1);
  }
  c.source = observed_pos;
  return c;
}
