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

#ifndef __HASH_H__
#define __HASH_H__

#include "klib/kvec.h" // C dynamic vector
#include "klib/khash.h" // C hash table/dictionary

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

typedef kvec_t(uint8_t) byteVec;


// creates uint32(qgram hash):kvec<readNum,pos> hash
KHASH_MAP_INIT_INT(qgramHash, matchVec);

// creates uint32(target read id):kvec<qpos,tpos> hash
KHASH_MAP_INIT_INT(matchHash, pairVec);

// similar to klib's khash.h, except we explicitly limit the length - it doesn't have to be null-terminated
static kh_inline khint_t qgram_hash(uint8_t *s, int k)
{
  int i;
  khint_t h = (khint_t)*s;
  if (h) for (i = 1 ; i < k; i++) h = (h << 5) - h + (khint_t)*(s+i);
  return h;
}

int hash_rmap(rmap map, int q, int threshold, int max_qgrams, int readLimit);

#endif /* __HASH_H__ */
