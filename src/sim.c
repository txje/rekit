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

#define PI 3.14159265358979323


// from https://en.wikipedia.org/wiki/Box-Muller_transform
float normal(float mu, float sigma) {
  const float two_pi = 2.0 * PI;
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

// inverse of the Cauchy CDF
float cauchy(float location, float scale) {
  int r = rand();
  /*
  fprintf(stderr, "rand: %d\n", r);
  fprintf(stderr, "randmax: %d\n", RAND_MAX);
  fprintf(stderr, "frac: %f\n", (r * (1.0/RAND_MAX)));
  fprintf(stderr, "tan(...): %f\n", tan(PI * (r * (1.0 / RAND_MAX) - 0.5)));
  */
  return scale * tan(PI * (r * (1.0 / RAND_MAX) - 0.5)) + location;
}

// end_idx is *not included* itself
u32Vec* bn_map(u32Vec *positions, int start_idx, int end_idx, uint64_t start_pos, uint32_t frag_len, float fn_rate, float fp_rate, float err_mean, float err_std, uint32_t resolution_min) {
  int j, k;

  int fp = round((end_idx - start_idx) * normal(fp_rate, 0.01));
  u32Vec fp_pos;
  kv_init(fp_pos);
  for(j = 0; j < fp; j++) {
    kv_push(uint32_t, fp_pos, (uint32_t)round((double)rand() / (double)RAND_MAX * frag_len));
  }
  ks_mergesort(uint32_t, kv_size(fp_pos), fp_pos.a, 0); // sort
  k = 0; // index into fp_pos

  // compute per-molecule uniform stretch by observed (query given ref) size / ref
  float uniform_stretch = (3014.8 + 0.955764 * frag_len) * normal(1.03025, 0.03273) / frag_len;
  //fprintf(stderr, "\nuniform stretch factor: %f\n", uniform_stretch);

  // now perform modifications for FN, FP, sizing error, and limited resolution
  u32Vec* modpos = (u32Vec*)malloc(sizeof(u32Vec));
  kv_init(*modpos);
  uint32_t last = 0;
  uint32_t last_stretched = 0;
  for(j = start_idx; j <= end_idx; j++) {
    uint32_t val;
    if(j < end_idx) {
      if(k < kv_size(fp_pos) && kv_A(fp_pos, k) < kv_A(*positions, j) - start_pos) {
        val = kv_A(fp_pos, k);
        k++;
        j--;
      } else {
        val = kv_A(*positions, j) - start_pos;
      }
    } else { // add a label for the end of the fragment (which DOES NOT correspond to a label site)
      val = frag_len;
    }
    // apply Cauchy-distributed inter-label error
    float c = cauchy(err_mean, err_std);
    while(c < 0) c = cauchy(err_mean, err_std); // under some error parameters, a proper cauchy random variable will end up with negative values, which we can't allow
    uint32_t f = last_stretched + (val - last) * uniform_stretch * c;

    // then include only fragments that exceed some minimum size (typically, ~1kb for Bionano)
    // and fall above FN rate
    if(((double)rand() / (double)RAND_MAX) > fn_rate || j == end_idx) { // this last position is the end of the molecule and can't be FN
      if(kv_size(*modpos) == 0 || f - last_stretched >= resolution_min) {
        kv_push(uint32_t, *modpos, f);
      } else { // if this label is too close to the last, use only the midpoint of the two
        kv_A(*modpos, kv_size(*modpos)-1) = last_stretched + (f - last_stretched) / 2;
      }
    }

    last = val;
    last_stretched = f;
  }

  return modpos;
}


cmap simulate_bnx(char* ref_fasta, char** motifs, size_t n_motifs, float frag_prob, float fn, float fp, float err_mean, float err_std, uint32_t resolution_min, float coverage) {

  float bimera_prob = 0.01;
  float trimera_prob = 0.0001;
  float quadramera_prob = 0.000001;

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

  fragVec ref_labels;
  kv_init(ref_labels);

  u32Vec ref_lens;
  kv_init(ref_lens);

  cstrVec ref_names;
  kv_init(ref_names);

  f = gzopen(ref_fasta, "r");
  seq = kseq_init(f);
  fprintf(stderr, "Loading FASTA file: %s\n", ref_fasta);

  int l;
  uint32_t ref = 0; // reference seq ID
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //fprintf(stderr, "Reading sequence '%s' (%i bp).\n", seq->name.s, l);
    genome_size = genome_size + seq->seq.l;

    int nlimit = 100; // break strings of Ns at least 100bp long

    // do digestion + shearing with some probability
    u32Vec *positions = (u32Vec*)malloc(sizeof(u32Vec));
    kv_init(*positions);
    int res = digest(seq->seq.s, l, motifs, n_motifs, 1, 0, nlimit, positions);
    //fprintf(stderr, "  has %u labels.\n", kv_size(*positions));

    char* n = malloc((seq->name.l+1)*sizeof(char));
    n[seq->name.l] = '\0';
    strncpy(n, seq->name.s, seq->name.l);
    kv_push(char*, ref_names, n);
    kv_push(u32Vec*, ref_labels, positions);
    kv_push(uint32_t, ref_lens, seq->seq.l); // I think this comes as a size_t, hopefully it will cast quietly...

    ref++;
  }
  kseq_destroy(seq);
  gzclose(f);

  uint64_t target_coverage = (uint64_t)((double)coverage * genome_size);
  uint64_t tot_covg = 0;
  fprintf(stderr, "Target bp: %llu (%fx coverage of %u bp genome)\n", target_coverage, coverage, genome_size);
  double chimera_prob;
  int chimera_parts = 0;
  u32Vec* prev_f;
  uint64_t bigrand;
  uint32_t ref_id;
  uint64_t pos;
  uint32_t frag_len;
  uint32_t last;
  uint32_t tmp;

  for(tot_covg = 0; tot_covg < target_coverage; ) {
    bigrand = rand(); // proof ourselves against 32-bit systems where our rand() will not necessarily be big enough 
    bigrand = (bigrand << 32) | rand();
    pos = bigrand % genome_size;
    //fprintf(stderr, "raw pos %lu\n", pos);
    for(i = 0; i < ref; i++) {
      if(pos < kv_A(ref_lens, i)) {
        ref_id = i;
        break;
      } else {
        pos = pos - kv_A(ref_lens, i);
        //fprintf(stderr, "skipping over ref %u of size %u\n", i, kv_A(ref_lens, i));
      }
    }
    //fprintf(stderr, "ref %u, pos %lu\n", ref_id, pos);
    ref_pos rp = {ref_id, (uint32_t)pos};
    frag_len = log(1 - (double)rand()/(double)RAND_MAX) / log(1 - frag_prob);
    //fprintf(stderr, "fragment length: %u\n", frag_len);
    if(pos + frag_len > kv_A(ref_lens, ref_id))
      frag_len = kv_A(ref_lens, ref_id) - pos;

    // find fragment start and end index in ref_pos
    for(i = 0; i < kv_size(*kv_A(ref_labels, ref_id)) && kv_A(*kv_A(ref_labels, ref_id), i) < pos; i++);
    for(j = i; j < kv_size(*kv_A(ref_labels, ref_id)) && kv_A(*kv_A(ref_labels, ref_id), j) < pos+frag_len; j++);
    //fprintf(stderr, "labels %d -> %d\n", i, j);

    u32Vec *f = bn_map(kv_A(ref_labels, ref_id), i, j, pos, frag_len, fn, fp, err_mean, err_std, resolution_min);
    //fprintf(stderr, "mapping done, got %u fragments\n", kv_size(*f));

    // reverse fragments randomly to represent opposite strand (labels are already strand-agnostic)
    if(rand()%2 == 0) {
      uint32_t len = kv_A(*f, kv_size(*f)-1);
      for(i = 0; i < (kv_size(*f)-1)/2; i++) { // we leave the last label alone since it represents the molecule length
        tmp = kv_A(*f, i);
        kv_A(*f, i) = len - kv_A(*f, kv_size(*f)-2-i); // the rest get reversed in order and adjusted to remain monotonically increasing
        kv_A(*f, kv_size(*f)-2-i) = len - tmp;
      }
    }

    // apply chimerism
    if(chimera_parts == 0) {
      //fprintf(stderr, "new chimera check\n");
      chimera_prob = (double)rand() / (double)RAND_MAX;
      chimera_parts = (chimera_prob < quadramera_prob ? 4 : (chimera_prob < trimera_prob ? 3 : (chimera_prob < bimera_prob ? 2 : 1)));
      //fprintf(stderr, "  chimera parts: %d\n", chimera_parts);

      kv_push(ref_pos, frag_positions, rp); // does not record the other parts (if any) of a chimera
      if(chimera_parts == 1) { // normal, never chimeric
        //fprintf(stderr, "  pushing singleton\n");
        kv_push(u32Vec*, fragments, f);
        chimera_parts = 0;
      } else {
        //fprintf(stderr, "  beginning of new chimera\n");
        prev_f = f;
        chimera_parts--;
      }
    } else {
      //fprintf(stderr, "  adding part %d to chimera\n", chimera_parts);
      // shift all of the f1 positions over to append to f0
      last = kv_A(*prev_f, kv_size(*prev_f)-1); // largest value
      for(j = 0; j < kv_size(*f); j++) {
        kv_push(uint32_t, *prev_f, kv_A(*f, j) + last);
      }

      kv_destroy(*f);
      free(f);
      if(chimera_parts == 1)
        kv_push(u32Vec*, fragments, prev_f);
      chimera_parts--;
    }

    tot_covg = tot_covg + frag_len;
    //fprintf(stderr, "-- running total: %u of %u\n", tot_covg, target_coverage);
  }

  // if last was an incomplete chimera, just end it and add it
  if(chimera_parts > 0)
    kv_push(u32Vec*, fragments, prev_f);

  cmap c;
  init_cmap(&c);
  c.n_rec_seqs = n_motifs;
  c.rec_seqs = motifs;
  fprintf(stderr, "Sampled %u fragments\n", kv_size(fragments));
  for(i = 0; i < kv_size(fragments); i++) {
    //fprintf(stderr, "adding map %d of size %u\n", i, kv_size(*kv_A(fragments, i)));
    add_map(&c, i+1, kv_A(fragments, i)->a, kv_size(*kv_A(fragments, i)), 1);
  }
  c.source = frag_positions;
  return c;
}
