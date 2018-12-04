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
  ref -->
t +----------------------
r |   ----> DEL
a |  |
c |  |   \
e |  |    \
  |  v     \
| |  INS  [MIS]MATCH
| |
v |
  |
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "klib/kvec.h" // C dynamic vector
#include "dtw.h"
#include <immintrin.h> // x86 intrinsics (SSE/AVX vector operations)

#define aln_gt(a,b) ((a).score > (b).score)
KSORT_INIT(aln_cmp, result, aln_gt)

/*
 * Overlap dynamic programming (time warping) alignment
 *
 * First row and column are initialized to zero, and alignment must reach either the last row or column
 */
result dtw(uint32_t* query, uint32_t* target, size_t qlen, size_t tlen, int8_t ins_score, int8_t del_score, float neutral_deviation, int rev) {
  result res;
  res.qrev = rev;

  if(tlen == 0 || qlen == 0) {
    res.failed = 1;
    res.score = -1;
    // other fields are unset and unreliable
    return res;
  }

  // must be dynamically allocated to prevent stack overflows (on *some* systems)
  // to make this more efficient for our uses, they could be allocated once for the largest sequences and reused
  // TODO: score and size matrices actually only need to keep 2 rows
  float* score_matrix[qlen+1];
  uint8_t* direction_matrix[qlen+1];
  size_t* q_cum_size[qlen+1]; // keeps track of accumulated sizes of indel'd fragments
  size_t* t_cum_size[qlen+1]; // keeps track of accumulated sizes of indel'd fragments
  int i;
  for(i = 0; i <= qlen; i++) {
    score_matrix[i] = (float*)malloc((tlen+1) * sizeof(float));
    direction_matrix[i] = (uint8_t*)malloc((tlen+1) * sizeof(uint8_t));
    q_cum_size[i] = (size_t*)malloc((tlen+1) * sizeof(size_t));
    t_cum_size[i] = (size_t*)malloc((tlen+1) * sizeof(size_t));
  }

  int x, y, qy; // qy will be y in fw direction, and qlen-1-y if rev

  // initialize first row and first column for local alignment
  for(y = 1; y <= qlen; y++) {
    score_matrix[y][0] = 0; // LOW for full-query, partial-target alignment
    q_cum_size[y][0] = 0;
    t_cum_size[y][0] = 0;
  }
  for(x = 0; x <= tlen; x++) {
    score_matrix[0][x] = 0;
    q_cum_size[0][x] = 0;
    t_cum_size[0][x] = 0;
  }

  float match, qmatch, tmatch, qtmatch, ins, del; // qmatch and tmatch are different kinds of matches were it accounts for only the extra q_size or t_size

  for(y = 0; y < qlen; y++) {
    qy = rev ? qlen-1-y : y;
    for(x = 0; x < tlen; x++) {
      // resetting any negative values to 0 is what makes this local alignment - if you don't do that it will be at least semi-global
      qtmatch = score_matrix[y][x] + score(q_cum_size[y][x] + query[qy], t_cum_size[y][x] + target[x], neutral_deviation) + 0.2;
      tmatch = score_matrix[y][x] + score(query[qy], t_cum_size[y][x] + target[x], neutral_deviation) + 0.1;
      qmatch = score_matrix[y][x] + score(q_cum_size[y][x] + query[qy], target[x], neutral_deviation) + 0.1;
      match = score_matrix[y][x] + score(query[qy], target[x], neutral_deviation);
      // basically, you get a bonus for using the leftover size from skipped fragments, but you don't have to
      match = match > qmatch && match > tmatch && match > qtmatch ? match : (qtmatch > qmatch && qtmatch > tmatch ? qtmatch : (qmatch > tmatch ? qmatch : tmatch));
      ins = score_matrix[y][x+1] + ins_score;
      del = score_matrix[y+1][x] + del_score;

      //printf("%i,%i: match %f, ins %f, del %f, score(%f,%f) %f\n", x, y, match, ins, del, query[qy], target[x], score(query[qy], target[x], neutral_deviation));

      // pick the highest-scoring move to make
      if(match >= ins && match >= del) {
        score_matrix[y+1][x+1] = match;
        direction_matrix[y+1][x+1] = MATCH;
        q_cum_size[y+1][x+1] = 0;
        t_cum_size[y+1][x+1] = 0;
      } else if(ins >= del) {
        score_matrix[y+1][x+1] = ins;
        direction_matrix[y+1][x+1] = INS;
        q_cum_size[y+1][x+1] = q_cum_size[y][x+1] + query[qy];
      } else {
        score_matrix[y+1][x+1] = del;
        direction_matrix[y+1][x+1] = DEL;
        t_cum_size[y+1][x+1] = t_cum_size[y+1][x] + target[x];
      }
    }
  }

  /*
  // print query
  printf("query: ");
  for(x = 0; x < 10; x++) {
    printf("\t%f", query[x]);
  }
  printf("\n");

  // print target
  printf("target: ");
  for(x = 0; x < 10; x++) {
    printf("\t%f", target[x]);
  }
  printf("\n");

  // print matrix
  for(y = 0; y < 10; y++) {
    printf("\n");
    for(x = 0; x < 10; x++) {
      printf("\t%f", score_matrix[y][x]);
    }
  }
  */

  // compute maximum score position (anywhere in last row or column to capture overlaps)
  int max_x = 0, max_y = 0;
  for(x = 1; x <= tlen; x++) {
    if(score_matrix[qlen][x] > score_matrix[max_y][max_x]) {
      max_x = x;
      max_y = qlen;
    }
  }
  for(y = 1; y <= qlen; y++) {
    if(score_matrix[y][tlen] > score_matrix[max_y][max_x]) {
      max_x = tlen;
      max_y = y;
    }
  }

  x = max_x;
  y = max_y;
  kv_init(res.path);
  while(y > 0 && x > 0) {
    kv_push(uint8_t, res.path, direction_matrix[y][x]);
    if(direction_matrix[y][x] == MATCH) {
      x--;
      y--;
    } else if(direction_matrix[y][x] == INS) {
      y--;
    } else if(direction_matrix[y][x] == DEL) {
      x--;
    }
  }

  res.score = score_matrix[max_y][max_x];
  res.qstart = y;
  res.qend = max_y;
  res.tstart = x;
  res.tend = max_x;
  res.failed = 0;
  // end positions are INCLUSIVE
  
  //printf("max_y: %i, max_x: %i\n", max_y, max_x);
  //printf("DTW score: %f\n", score_matrix[max_y][max_x]);

  // free these up in case we'll be doing this repeatedly
  for(i = 0; i <= qlen; i++) {
    free(score_matrix[i]);
    free(direction_matrix[i]);
    free(q_cum_size[i]);
    free(t_cum_size[i]);
  }

  return res;
}
