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

#ifndef __DTW_H__
#define __DTW_H__

#include "klib/kvec.h" // C dynamic vector
#include "rmap.h"

typedef kvec_t(uint8_t) pathvec;

typedef struct aln_result {
  float score;
  uint32_t qstart;
  uint32_t qend;
  uint32_t tstart;
  uint32_t tend;
  uint8_t failed; // boolean flag
} result;

static uint8_t MATCH = 0, INS = 1, DEL = 2;

static float LOW = -1e38; // the order of the lowest possible 4-byte (single-precision) float


static float score(float a, float b, float neutral_deviation) {
  float diff = (a > b ? (a - b) : (b - a));
  // score is:
  // 1 if fragments are the same size
  // 0 if fragments differ by neutral deviation
  // -1 if size difference is 2x neutral deviation
  return 1 - (diff / neutral_deviation);
}

int dtw_rmap(rmap map, int threshold);
result align(float* query, float* target, size_t qlen, size_t tlen, pathvec *path, int8_t ins_score, int8_t del_score, float neutral_deviation);

#endif /* __DTW_H__ */
