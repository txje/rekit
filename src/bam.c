#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "htslib/sam.h"
#include "cmap.h"

/*
 * Jeremy Wang
 * 20180817
*/

cmap get_cmap_from_bam(char* bam_file, int covg_threshold) {

  cmap c;
  init_cmap(&c);

  // load BAM file
  //
  samFile *bam;
  bam_hdr_t *header;
  bam1_t *aln;
  int ret_val;

  if(strcmp(bam_file,  "-") == 0) {
    bam = sam_open("-", "r");
    //fprintf(stderr, "Reading from stdin...\n");
  } else {
    bam = sam_open(bam_file, "rb");
    //fprintf(stderr, "Reading from bam file '%s'...\n", bam_file);
  }
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", bam_file);
    return c;
  }
  header = sam_hdr_read(bam);
  if (header == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_file);
    return c;
  }
  // construct array from reference information so that we can look it up with read.tid
  int **covg = malloc(sizeof(int*) * header->n_targets);
  int *rlen_array = malloc(sizeof(int) * header->n_targets); // has to be a plain int because that's what kseq gives out

  int i, j;
  for (i = 0; i < header->n_targets; i++) {
    covg[i] = calloc(header->target_len[i] / 100 + 1, sizeof(int));

    // sequence length
    rlen_array[i] = header->target_len[i];
  }

  aln = bam_init1();

  int32_t read_len = 0;

  uint32_t ct = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {

    if (aln->core.flag & 4+256+2048) { // unmapped, secondary, or supplementary alignment
      continue;
    }

    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;

    covg[tid][pos/100]++;
    if(pos/100 != endpos/100) {
      covg[tid][endpos/100]++;
    }

    ct++;
    if(ct % 1000000 == 0) {
      //fprintf(stderr, "%d reads processed\n", ct);
    }
  }

  uint64_t bin;
  // go through seqs, accumulate per-sample stats...
  //fprintf(stderr, "chrom\tpos\tcoverage\n");
  uint64_t tot_bin_covg, tot_bin_pos, last_bin;
  u32Vec pos;
  for (i = 0; i < header->n_targets; i++) {
    kv_init(pos);

    tot_bin_covg = 0;
    tot_bin_pos = 0;
    last_bin = 0;
    for(bin = 0; bin < rlen_array[i]/100+1; bin++) {
      if(covg[i][bin] > 0) {
        if(bin == last_bin + 1) {
          tot_bin_covg += covg[i][bin];
          tot_bin_pos += (bin*100 * covg[i][bin]);
        } else{
          if (tot_bin_covg >= covg_threshold) {
            //fprintf(stderr, "%s\t%d\t%d\n", header->target_name[i], tot_bin_pos/tot_bin_covg+50, tot_bin_covg);
            kv_push(uint32_t, pos, tot_bin_pos/tot_bin_covg+50);
          }
          tot_bin_covg = covg[i][bin];
          tot_bin_pos = (bin*100 * covg[i][bin]);
        }
        last_bin = bin;
      }
    }
    if(tot_bin_covg > covg_threshold) { // last bin, if necessary
      //fprintf(stderr, "%s\t%d\t%d\n", header->target_name[i], tot_bin_pos/tot_bin_covg+50, tot_bin_covg);
      kv_push(uint32_t, pos, tot_bin_pos/tot_bin_covg+50);
    }
    // put in the end of the chromosome
    kv_push(uint32_t, pos, rlen_array[i]);

    add_map(&c, i+1, pos.a, kv_size(pos), 1);
    kv_destroy(pos); // the values were copied to an array of labels, so we can free these positions
  }

  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return c;
  }

  bam_destroy1(aln);

  return c;
}
