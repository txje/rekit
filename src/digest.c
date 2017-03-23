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
#include "klib/kvec.h"
#include "digest.h"

int digest(char *seq, size_t seq_len, char **motifs, size_t n_motifs, float digest_rate, float shear_rate, int nlimit, u32Vec *sizes) {
  srand(time(NULL));

  // find first non-N character in seq
  int start_pos = 0;
  while(start_pos < seq_len && seq[start_pos] == 'N') {
    start_pos++;
  }
  int last = start_pos;
  //printf("Non-N start position: %d\n", start_pos);

  if(start_pos >= seq_len) {
    return 1; // all N sequence
  }

  // go through sequence
  int nstart = -1;
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

  for(i = start_pos; i < seq_len; i++) {
    //printf("At pos %d\n", i);

    // cut at strings of Ns at least <nlimit> long
    // and EXCISE THE N PORTION
    if(seq[i] == 'N') {
      if(nstart == -1) {
        nstart = i;
      }
      continue;
    } else {
      if(nstart != -1 && i - nstart >= nlimit) {
        if(nstart - last > 0) {
          kv_push(uint32_t, *sizes, (uint32_t)(nstart - last));
        }
        last = i;
      }
      nstart = -1;
    }

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
            printf("Error: Invalid character encountered at seq pos %d: %c\n", i+j, seq[i+j]);
            return 2;
        }
      }
      //printf("Comparing %s and %s (rc) to motif %s (%dbp)\n", mers[0][m], mers[1][m], motifs[m], mlen);

      float digest_chance = (double)rand() / (double)RAND_MAX;
      float shear_chance = (double)rand() / (double)RAND_MAX;
      //printf("Random number %f < rate %f?\n", rnd, digest_rate);
      if(i - last > 0 && (shear_chance < shear_rate || ((strcmp(motifs[m], mers[0][m]) == 0 || strcmp(motifs[m], mers[1][m]) == 0) && digest_chance < digest_rate))) {
        kv_push(uint32_t, *sizes, (uint32_t)(i - last));
        last = i;
        break; // stop checking motifs
      }
    }
  }

  if((nstart == -1 || i - nstart < nlimit) && i - last > 0) {
    kv_push(uint32_t, *sizes, (uint32_t)(i - last));
  } else if(nstart - last > 0){
    kv_push(uint32_t, *sizes, (uint32_t)(nstart - last));
  }

  // ended without incident
  return 0;
}
