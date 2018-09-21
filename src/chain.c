#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "klib/ksort.h"
#include "hash.h"

typedef struct score_pos {
  int score;
  int pos;
  int prev;
  uint8_t used;
} score_pos;

#define score_pos_gt(a,b) ((a).score > (b).score)
KSORT_INIT(score_pos_cmp, score_pos, score_pos_gt)

// input anchors should be sorted by tpos increasing
pairVec* chain(pairVec* anchors, int max_chains, int match_score, int max_gap, int min_chain_length) {

  score_pos* scores = malloc(sizeof(score_pos) * kv_size(*anchors));
  scores[0].score = match_score;
  scores[0].pos = 0;
  scores[0].prev = -1;
  scores[0].used = 0;

  score_pos max_score;
  int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
  int h = 50; // number of previous anchors to check
  for(i = 1; i < kv_size(*anchors); i++) {
    //fprintf(stderr, "computing score at pos %d\n", i);

    max_score.score = match_score;
    max_score.pos = i;
    max_score.prev = -1; // no predecessor
    max_score.used = 0;

    for(j = i < h ? 0 : i-h; j < i; j++) {
      qdiff = (int)kv_A(*anchors, i).qpos - (int)kv_A(*anchors, j).qpos;
      tdiff = (int)kv_A(*anchors, i).tpos - (int)kv_A(*anchors, j).tpos;
      if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) { // we often have multiple hits to the same target pos, so we can't let them chain with either other
        continue;
      }
      diffdiff = (qdiff > tdiff ? qdiff - tdiff : tdiff - qdiff);
      gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff));
      qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
      score = scores[j].score + (qdiff > match_score ? match_score : qdiff) - gap_cost;
      if(score > max_score.score) {
        max_score.score = score;
        max_score.prev = j;
      }
    }
    scores[i] = max_score;
    //fprintf(stderr, "  setting max score to %d from pos %d\n", max_score.score, max_score.prev);
  }

  // copy unsorted scores:
  score_pos* anchor_scores = malloc(sizeof(score_pos) * kv_size(*anchors));
  memcpy(anchor_scores, scores, sizeof(score_pos) * kv_size(*anchors));
  //sort scores decreasing
  ks_mergesort(score_pos_cmp, kv_size(*anchors), scores, 0);

  // build non-overlapping chains from highest to lowest score
  // and merge overlapping chains
  // we're going to use these vectors in a way they're not really supposed to be used (like a normal array) so that this function will be backwards-compatible with the pairVec
  pairVec* chains = malloc(sizeof(pairVec) * max_chains);
  int c; // chain index
  i = 0; // scores index
  for(c = 0; c < max_chains && i < kv_size(*anchors); c++) {
    //fprintf(stderr, "building chain %d\n", c);
    int chain_len = 0;
    int chain_pos = scores[i].pos; // get the pos which should be index into the anchor_scores array
    while(anchor_scores[chain_pos].used == 0) {
      //fprintf(stderr, "backtracking from %d to %d\n", chain_pos, scores[chain_pos].prev);
      chain_len++;
      if(anchor_scores[chain_pos].prev == -1)
        break;
      chain_pos = anchor_scores[chain_pos].prev;
    }
    if(chain_len < min_chain_length) {
      c--;
      i++;
      continue;
    }
    kv_init(chains[c]);
    kv_resize(posPair, chains[c], chain_len);
    chains[c].n = chain_len; // set the size explicitly, then we'll set the values explicitly
    //fprintf(stderr, "creating chain %d of length %d\n", c, kv_size(chains[c]));
    chain_pos = scores[i].pos;
    for(j = 0; j < chain_len; j++) {
      anchor_scores[chain_pos].used = 1;
      kv_A(chains[c], chain_len-1-j) = kv_A(*anchors, anchor_scores[chain_pos].pos);
      chain_pos = anchor_scores[chain_pos].prev;
    }
    i++;
  }
  if(c < max_chains) kv_init(chains[c]); // an empty vector will indicate the end of the chains array if there are fewer than max_chains results

  return chains;
}

