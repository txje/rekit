#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "klib/ksort.h"
#include "klib/khash.h"
#include "klib/kvec.h"
#include "hash.h"
#include "chain.h"

#define score_pos_gt(a,b) ((a).score > (b).score)
KSORT_INIT(score_pos_cmp, score_pos, score_pos_gt)

#define pos_pair_lt(a,b) ((a).tpos < (b).tpos)
KSORT_INIT(pos_pair_cmp, posPair, pos_pair_lt)

// input anchors should be sorted by tpos increasing
chain* do_chain(khash_t(matchHash) *hits, int max_chains, int match_score, int max_gap, int min_chain_length) {

  // assess hits for each target
  pairVec anchors;
  khint_t bin;
  uint32_t target;

  scoreVec scores;
  kv_init(scores);
  score_pos s;

  int i, score, qdiff, tdiff, diffdiff, gap_cost, j, best_j;
  int h = 50; // number of previous anchors to check
  int ref_offset = 0; // offset into the scores vector of the current target - sets the backward limit for finding chained anchors

  // iterate through hits for each target, and append them to the same scores vector
  for (bin = kh_begin(hits); bin != kh_end(hits); ++bin) {
    if (!kh_exist(hits, bin)) continue;
    target = kh_key(hits, bin);
    anchors = kh_val(hits, bin);

    //sort anchor pairs by target pos increasing
    ks_mergesort(pos_pair_cmp, kv_size(anchors), anchors.a, 0);

    s.score = match_score;
    s.pos = ref_offset;
    s.prev = -1;
    s.used = 0;
    s.ref = target;
    kv_push(score_pos, scores, s); // this actually copies the score_pos struct values, so we can reuse 's'

    for(i = 1; i < kv_size(anchors); i++) {

      s.score = match_score;
      s.pos = i;
      s.prev = -1; // no predecessor

      for(j = i < h ? 0 : i-h; j < i; j++) {
        qdiff = (int)kv_A(anchors, i).qpos - (int)kv_A(anchors, j).qpos;
        tdiff = (int)kv_A(anchors, i).tpos - (int)kv_A(anchors, j).tpos;
        if(tdiff <= 0 || qdiff <= 0 || qdiff > max_gap || tdiff > max_gap) { // we often have multiple hits to the same target pos, so we can't let them chain with either other
          continue;
        }
        diffdiff = (qdiff > tdiff ? qdiff - tdiff : tdiff - qdiff);
        gap_cost = (diffdiff == 0 ? 0 : 0.01 * match_score * diffdiff + 0.5 * log2(diffdiff)); // affine gap penalty a la minimap2
        qdiff = (qdiff > tdiff ? tdiff : qdiff); // now not qdiff but the minimum difference
        score = kv_A(scores, j + ref_offset).score + (qdiff > match_score ? match_score : qdiff) - gap_cost;
        if(score > s.score) {
          s.score = score;
          s.prev = j + ref_offset;
        }
      }
      kv_push(score_pos, scores, s);
    }
    ref_offset = kv_size(scores);
  }

  // copy unsorted scores:
  score_pos* anchor_scores = malloc(sizeof(score_pos) * kv_size(scores));
  memcpy(anchor_scores, scores.a, sizeof(score_pos) * kv_size(scores));
  //sort scores decreasing
  ks_mergesort(score_pos_cmp, kv_size(scores), scores.a, 0);

  // build non-overlapping chains from highest to lowest score
  chain* chains = malloc(sizeof(chain) * max_chains); // * 1.05);
  int c; // chain index
  i = 0; // scores index
  for(c = 0; c < max_chains && i < kv_size(scores); c++) {
    //fprintf(stderr, "building chain %d\n", c);
    int chain_len = 0;
    int chain_pos = kv_A(scores, i).pos; // get the pos which should be index into the anchor_scores array
    // get anchors associated with this score/ref
    bin = kh_get(matchHash, hits, kv_A(scores, i).ref);
    if(bin == kh_end(hits)) { // key not found, *shouldn't* happen
      fprintf(stderr, "something went very wrong with the chain computation!");
      return (chain*)NULL;
    }
    anchors = kh_val(hits, bin);
    while(anchor_scores[chain_pos].used == 0) {
      //printf("backtracking from %d (%d:%d) to %d\n", chain_pos, kv_A(anchors, anchor_scores[chain_pos].pos).qpos, kv_A(anchors, anchor_scores[chain_pos].pos).tpos, anchor_scores[chain_pos].prev);
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
    kv_init(chains[c].anchors);
    chains[c].ref = kv_A(scores, i).ref;
    //printf("chain %d (score %d, %d anchors) is from ref %u anchor %d\n", c, kv_A(scores, i).score, chain_len, chains[c].ref, kv_A(scores, i).pos);
    kv_resize(posPair, chains[c].anchors, chain_len);
    chains[c].anchors.n = chain_len; // set the size explicitly, then we'll set the values explicitly
    //fprintf(stderr, "creating chain %d of length %d\n", c, kv_size(chains[c]));
    chain_pos = kv_A(scores, i).pos;
    for(j = 0; j < chain_len; j++) {
      anchor_scores[chain_pos].used = 1;
      kv_A(chains[c].anchors, chain_len-1-j) = kv_A(anchors, anchor_scores[chain_pos].pos);
      chain_pos = anchor_scores[chain_pos].prev;
    }
    i++;
  }
  if(c < max_chains) kv_init(chains[c].anchors); // an empty vector will indicate the end of the chains array if there are fewer than max_chains results

  free(anchor_scores);
  kv_destroy(scores);
  return chains;
}

