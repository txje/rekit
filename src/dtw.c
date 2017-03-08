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
#include "klib/kvec.h" // C dynamic vector
#include "dtw.h"

result align(float* query, float* target, uint32_t qlen, uint32_t tlen, pathvec *path, int8_t ins_score, int8_t del_score, float neutral_deviation) {

  if(tlen == 0 || qlen == 0) {
    result res;
    res.failed = 1;
    // other fields are unset and unreliable
    return res;
  }

  // must be dynamically allocated to prevent stack overflows (on *some* systems)
  float* score_matrix[qlen+1];
  uint8_t* direction_matrix[qlen+1];
  int i;
  for(i = 0; i <= qlen; i++) {
    score_matrix[i] = (float*)malloc((tlen+1) * sizeof(float));
    direction_matrix[i] = (uint8_t*)malloc((tlen+1) * sizeof(uint8_t));
  }

  int x, y;

  // initialize matrix for full-query, partial-target alignment
  for(y = 1; y <= qlen; y++) {
    score_matrix[y][0] = LOW;
  }
  for(x = 0; x <= tlen; x++) {
    score_matrix[0][x] = 0;
  }

  float match, ins, del;

  for(y = 0; y < qlen; y++) {
    for(x = 0; x < tlen; x++) {
      // match
      match = score_matrix[y][x] + score(query[y], target[x], neutral_deviation);
      // ins
      ins = score_matrix[y][x+1] + ins_score;
      // del
      del = score_matrix[y+1][x] + del_score;
      //printf("%i,%i: match %f, ins %f, del %f, score(%f,%f) %f\n", x, y, match, ins, del, query[y], target[x], score(query[y], target[x], neutral_deviation));

      // compare
      if(match >= ins && match >= del) {
        score_matrix[y+1][x+1] = match;
        direction_matrix[y+1][x+1] = MATCH;
      } else if(ins >= del) {
        score_matrix[y+1][x+1] = ins;
        direction_matrix[y+1][x+1] = INS;
      } else {
        score_matrix[y+1][x+1] = del;
        direction_matrix[y+1][x+1] = DEL;
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

  // compute maximum score position (anywhere in last row)
  int max_x = 0, max_y = 0;
  max_x = 0;
  max_y = qlen;
  for(x = 1; x <= tlen; x++) {
    if(score_matrix[max_y][x] > score_matrix[max_y][max_x]) {
      max_x = x;
    }
  }

  x = max_x;
  y = max_y;
  while(y > 0) {
    kv_push(uint8_t, *path, direction_matrix[y][x]);
    if(direction_matrix[y][x] == MATCH) {
      x--;
      y--;
    } else if(direction_matrix[y][x] == INS) {
      y--;
    } else if(direction_matrix[y][x] == DEL) {
      x--;
    }
  }

  result res;
  res.score = score_matrix[max_y][max_x];
  res.qstart = y;
  res.qend = max_y;
  res.tstart = x;
  res.tend = max_x;
  // end positions are INCLUSIVE
  
  //printf("max_y: %i, max_x: %i\n", max_y, max_x);
  //printf("DTW score: %f\n", score_matrix[max_y][max_x]);

  // free these up in case we'll be doing this repeatedly
  for(i = 0; i <= qlen; i++) {
    free(score_matrix[i]);
    free(direction_matrix[i]);
  }

  return res;
}
