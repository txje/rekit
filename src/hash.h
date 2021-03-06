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

#include <math.h>
#include "klib/kvec.h" // C dynamic vector
#include "klib/khash.h" // C hash table/dictionary
#include "klib/ksort.h"
#include "cmap.h"

#ifndef __HASH_H__
#define __HASH_H__

// MAXIMUM # READS = 2^31 (~2bn)
typedef struct {
  uint32_t readNum; // last bit 0 if fw, 1 if rv (same as alternating fw/rv) [fw read n = n*2, rv read n = n*2+1]
  // in practice, we only hash the fw strand, and do lookup by fw and rv, so right now these will always end in 0 (readNum%2 == 0)
  uint32_t pos;
} readPos;

typedef struct {
  uint32_t qpos;
  uint32_t tpos;
} posPair;

typedef kvec_t(readPos) matchVec;

typedef kvec_t(posPair) pairVec;

// creates uint32(qgram hash):kvec<readNum,pos> hash
KHASH_MAP_INIT_INT(qgramHash, matchVec);

// creates uint32(target read id):kvec<qpos,tpos> hash
KHASH_MAP_INIT_INT(matchHash, pairVec);

static kh_inline khint_t xratio_hash(label* labels, int bins, int skip) {
  //fprintf(stderr, "%d %d %d %d\n", labels[0].position, labels[1].position, labels[2].position, labels[3].position);
  float cr;
  if(skip == 1)
    cr = (((float)labels[3].position - labels[0].position) * ((float)labels[4].position - labels[2].position)) / (((float)labels[3].position - labels[2].position) * ((float)labels[4].position - labels[0].position)); // cross-ratio
  else if(skip == 2)
    cr = (((float)labels[3].position - labels[0].position) * ((float)labels[4].position - labels[1].position)) / (((float)labels[3].position - labels[1].position) * ((float)labels[4].position - labels[0].position)); // cross-ratio
  else if(skip == 3)
    cr = (((float)labels[2].position - labels[0].position) * ((float)labels[4].position - labels[1].position)) / (((float)labels[2].position - labels[1].position) * ((float)labels[4].position - labels[0].position)); // cross-ratio
  else if(skip == 4)
    cr = (((float)labels[2].position - labels[0].position) * ((float)labels[3].position - labels[1].position)) / (((float)labels[2].position - labels[1].position) * ((float)labels[3].position - labels[0].position)); // cross-ratio
  //fprintf(stderr, "x-ratio: %f\n", cr);
  // cross-ratio CDF (per https://hal.inria.fr/inria-00590012/document), with some modifications
  // since all of our points are in increasing order, they are actually modeled by only the F1 portion of the CDF
  float crcdf = (1.0/2 + (cr*(1-cr)*log((cr-1)/cr) - cr + 1.0/2)) * 2;
  khint_t h = (khint_t)(crcdf * bins + bins * (labels[skip < 4 ? 4 : 3].position - labels[0].position)/2000); // the number of bins used here will impact the accuracy signficantly
  return h;
}

// similar to klib's khash.h, except we explicitly limit the length - it doesn't have to be null-terminated
static kh_inline khint_t qgram_hash(uint8_t *s, int k, int skip, uint32_t l) {
  int i, j = 0;
  khint_t h = 0;
  for (i = 0; i < k; i++) {
    //if(i != skip) {
      // l is a bit vector representing whether to ceil (+1) each value
      // if skipping the next label, add the two adjacent sizes together
      //h = (h << 5) - h + (khint_t)*(s+i) + (i+1==skip ? (khint_t)*(s+i+1) : 0) + (l>>j++ & 1);
      h = (h << 5) - h + (khint_t)*(s+i) + (l>>j++ & 1);
    //}
  }
  return h;
}

int hash_cmap(cmap b, cmap c, FILE* o, int q, int chain_threshold, float dtw_threshold, int max_qgrams, int readLimit, int bin_size, int resolution_min, int min_labels, int start_mol, int end_mol);

uint32_t* u32_get_fragments(label* labels, size_t n_labels, int bin_size, int rev);

#endif /* __HASH_H__ */
