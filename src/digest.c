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

#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "klib/kvec.h"
#include "digest.h"
#include "cmap.h"

#ifndef _kseq_
#define _kseq_

#include "klib/kseq.h"

// init kseq struct
KSEQ_INIT(gzFile, gzread)

#endif

int digest(char *seq, size_t seq_len, char **motifs, size_t n_motifs, float digest_rate, float shear_rate, int nlimit, u32Vec *positions) {
  srand(time(NULL));

  // go through sequence
  int i, j, m;

  // create a crazy triple ref holder for kmer/rc as we loop through the sequence
  char ***mers = malloc(sizeof(char**) * 2); // fw and rc
  for(i = 0; i < 2; i++) {
    mers[i] = (char**) malloc(sizeof(char*) * n_motifs);
    for(m = 0; m < n_motifs; m++) {
      mers[i][m] = (char*) malloc(sizeof(char) * (strlen(motifs[m]) + 1));
      mers[i][m][strlen(motifs[m])] = '\0';
    }
  }

  for(i = 0; i < seq_len; i++) {

    // see if current locus matches any motifs
    for(m = 0; m < n_motifs; m++) {
      size_t mlen = strlen(motifs[m]);
      if(i+mlen >= seq_len) {
        // not enough sequence left to match this motif
        continue;
      }

      // get kmer and reverse-complement to compare
      for(j = 0; j < mlen; j++) {
        switch(seq[i+j]){
          case 'A':
          case 'a':
            mers[0][m][j] = 'A';
            mers[1][m][mlen-j-1] = 'T';
            break;
          case 'C':
          case 'c':
            mers[0][m][j] = 'C';
            mers[1][m][mlen-j-1] = 'G';
            break;
          case 'G':
          case 'g':
            mers[0][m][j] = 'G';
            mers[1][m][mlen-j-1] = 'C';
            break;
          case 'T':
          case 't':
            mers[0][m][j] = 'T';
            mers[1][m][mlen-j-1] = 'A';
            break;
          case 'N':
          case 'n':
            mers[0][m][j] = 'N';
            mers[1][m][mlen-j-1] = 'N';
            break;
          default:
            fprintf(stderr, "Warning: Invalid character encountered at seq pos %d: %c\n", i+j, seq[i+j]);
            mers[0][m][j] = 'N';
            mers[1][m][mlen-j-1] = 'N';
        }
      }
      //printf("Comparing %s and %s (rc) to motif %s (%dbp)\n", mers[0][m], mers[1][m], motifs[m], mlen);

      /*
      float digest_chance = (double)rand() / (double)RAND_MAX;
      float shear_chance = (double)rand() / (double)RAND_MAX;
      //printf("Random number %f < rate %f?\n", rnd, digest_rate);
      if(shear_chance < shear_rate || ((strcmp(motifs[m], mers[0][m]) == 0 || strcmp(motifs[m], mers[1][m]) == 0) && digest_chance <= digest_rate)) {
        kv_push(uint32_t, *positions, (uint32_t)i);
        break; // stop checking motifs
      }
      */
      // do a perfect digestion - ignores FP and FN
      if (strcmp(motifs[m], mers[0][m]) == 0 || strcmp(motifs[m], mers[1][m]) == 0) {
        kv_push(uint32_t, *positions, (uint32_t)i);
        break;
      }
    }
  }

  kv_push(uint32_t, *positions, (uint32_t)i);

  return 0;
}

cmap digest_fasta(char* fasta_file, char** motifs, size_t n_motifs) {
  cmap c;
  init_cmap(&c);
  c.rec_seqs = motifs;
  c.n_rec_seqs = n_motifs;

  gzFile gzfp;
  kseq_t *seq;
  int l;

  gzfp = gzopen(fasta_file, "r");
  if(!gzfp) {
    fprintf(stderr, "File '%s' not found\n", fasta_file);
    return c;
  }
  seq = kseq_init(gzfp);

  int refid = 1;
  u32Vec labels;
  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //fprintf(stderr, "Reading %s (%i bp).\n", seq->name.s, l);
    kv_init(labels);
    digest(seq->seq.s, l, motifs, n_motifs, 1.0, 0.0, 1000000000, &labels); // 1.0 true digest rate, 0.0 random shear rate (perfect)
    //fprintf(stderr, "Produced %d labels\n", kv_size(labels));
    add_map(&c, refid++, labels.a, kv_size(labels), 1); // channel 1
  }

  gzclose(gzfp);
  kseq_destroy(seq);

  return c;
}
