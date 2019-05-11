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

#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

#include "cmap.h"
#include "bnx.h"
#include "hash.h"
//#include "lsh.h"
#include "sim.h"
#include "digest.h"
#include "bam.h"
#include "dtw.h"

void usage() {
  printf("Usage: rekit [command] [options]\n");
  printf("Commands:\n");
  printf("  align:    align BNX molecules to reference CMAP\n");
  printf("  dtw:      DTW-only align BNX molecules to reference CMAP\n");
  printf("  simulate: simulate molecules\n");
  printf("  digest:   in silico digestion\n");
  printf("  label:    produce alignment-based reference CMAP\n");
  printf("Options:\n");
  printf("  align    -bc\n");
  printf("  dtw      -bc\n");
  printf("  simulate -frx --break-rate --fn --fp --min-frag --stretch-mean --stretch-std --source-output\n");
  printf("  digest   -fr\n");
  printf("  label    -a\n");
  printf("    -b: bnx: A single BNX file containing molecules\n");
  printf("    -c: cmap: A single CMAP file\n");
  printf("    -f: fasta: Reference sequence to simulate from\n");
  printf("    -a: bam: BAM alignment file\n");
  printf("    -r: cutseq: Recognition/label site sequence\n");
  printf("    -q: Size of q-gram/k-mer to hash (default: 4)\n");
  printf("    -h: Number of hash functions to apply\n");
  //printf("    -e: Seed to random number generator\n");
  printf("    -t: Minimum number of q-gram/cross-ratio anchors in a chain (default: 1)\n");
  printf("    -m: max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored\n");
  printf("    -d: DTW score threshold to report alignment (default: 5)\n");
  printf("    -x: Simulated molecule coverage\n");
  printf("  simulate options:\n");
  printf("    --break-rate: Probability of genome fragmentation per locus (default: 0.000005)\n");
  printf("    --fn: Probability of missed label at true restriction site (default: 0.09893)\n");
  printf("    --fp: Probability of false-positive label (default: 0.07558)\n");
  printf("    --stretch-mean: Fragment stretch mean (default: 0.991385)\n");
  printf("    --stretch-std: Fragment stretch standard deviation (default: 0.033733)\n");
  printf("    --min-frag: Minimum detectable fragment size (default: 500)\n");
  printf("    -s, --source-output: Output the reference positions of the simulated molecules to the given file\n");
  printf("  label options:\n");
  printf("    --coverage-threshold: Read coverage required (in ~300bp window) to call a label site (default: 10)\n");
  printf("  align options:\n");
  printf("    --min-labels: Minimum molecule labels to align\n");
  printf("    --start-mol: Molecule ID to start at (for multithreading)\n");
  printf("    --end-mol: Molecule ID to end at (inclusive)\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "break-rate",             required_argument, 0, 0 },
  { "fn",                     required_argument, 0, 0 },
  { "fp",                     required_argument, 0, 0 },
  { "min-frag",               required_argument, 0, 0 },
  { "stretch-mean",           required_argument, 0, 0 },
  { "stretch-std",            required_argument, 0, 0 },
  { "coverage-threshold",     required_argument, 0, 0 },
  { "help",                   no_argument,       0, 0 },
  { "source-output",          required_argument, 0, 0 },
  { "bin-size",               required_argument, 0, 0 },
  { "min-labels",             required_argument, 0, 0 },
  { "start-mol",              required_argument, 0, 0 },
  { "end-mol",                required_argument, 0, 0 },
  { 0, 0, 0, 0}
};

int main(int argc, char *argv[]) {
  srand(time(NULL));

  char* bnx_file = NULL; // .bnx file path/name
  char* fasta_file = NULL; // .fasta file path/name
  char* cmap_file = NULL; // .cmap file path/name (either reference or consensus)
  char* bam_file = NULL; // .bam file path/name (aligned)
  char* restriction_seq = NULL; // restriction enzyme or label recognition sequence (must also be reverse complemented if not symmetrical)
  char* source_outfile = NULL; // output file for the truth/source positions
  int q = 5; // q-gram size (set to 5 to make sure when we go to hash we have 5 to make sets of 4-mers with each missing)
  int h = 10; // number of hashes
  int verbose = 0;
  int chain_threshold = 1;
  float dtw_threshold = 5;
  int seed = 0; // made this up
  int max_qgrams = 2000000000; // made this up
  int bin_size = 100; // # bins that x-ratios will be spread across, or divisor for fragment size binning
  int read_limit = -1; // just for testing
  int min_labels = 11; // a parameter, and this works well in practice
  int start_mol = 0;
  int end_mol = -1;

  float coverage = 0.0;
  int covg_threshold = 10;
  float break_rate = 0.000005; // one every 200Kb
  float fn = 0.09893; // based on empirical data for Saphyr DLE1
  float fp = 0.07558; // based on empirical data
  float min_frag = 500; // minimum reported gap between labels
  float stretch_mean = 0.991385; // these are based empirically on Cauchy distribution of NA12878 DLE1 data
  float stretch_std = 0.033733;

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "b:c:q:hf:r:t:m:vx:a:s:d:", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'b':
        bnx_file = optarg;
        break;
      case 'c':
        cmap_file = optarg;
        break;
      case 'q':
        q = atoi(optarg);
        break;
      case 'h':
        usage();
        return 0;
        break;
      case 'f':
        fasta_file = optarg;
        break;
      case 'r':
        restriction_seq = optarg;
        break;
      case 't':
        chain_threshold = atoi(optarg);
        break;
      case 'd':
        dtw_threshold = atof(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'x':
        coverage = atof(optarg);
        break;
      case 'a':
        bam_file = optarg;
        break;
      case 's':
        source_outfile = optarg;
        break;
      case '?':
        if (optopt == 'b' || optopt == 'c' || optopt == 'q' || optopt == 'r' || optopt == 'f' || optopt == 't' || optopt == 'm' || optopt == 'x' || optopt == 'a' || optopt == 's')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        //fprintf(stderr, "option %s\n", long_options[long_idx].name);
        if (long_idx == 0) break_rate = atof(optarg); // --break-rate
        else if (long_idx == 1) fn = atof(optarg); // --fn
        else if (long_idx == 2) fp = atof(optarg); // --fp
        else if (long_idx == 3) min_frag = atoi(optarg); // --min-frag
        else if (long_idx == 4) stretch_mean = atof(optarg); // --stretch-mean
        else if (long_idx == 5) stretch_std = atof(optarg); // --stretch-std
        else if (long_idx == 6) covg_threshold = atoi(optarg); // --coverage-threshold
        else if (long_idx == 7) {usage(); return 0;} // --help
        else if (long_idx == 8) source_outfile = optarg; // --source-output
        else if (long_idx == 9) bin_size = atoi(optarg); // --bin-size
        else if (long_idx == 10) min_labels = atoi(optarg); // --min-labels
        else if (long_idx == 11) start_mol = atoi(optarg)-1; // --start-mol, decrement to make it match 0-based indices instead of 1-based in BNX
        else if (long_idx == 12) end_mol = atoi(optarg)-1; // --end-mol
        break;
      default:
        usage();
        return 1;
    }
  }

  // CMAP test
  //write_cmap(&c, fopen("test.cmap", "w"));

  int index;
  char* command = NULL;
  for (index = optind; index < argc; index++) {
    if(index == optind) {
      command = argv[index];
    }
  }
  if(command == NULL) {
    usage();
    return 1;
  }

  cmap c;
  size_t n_frags;

  int i = 0;
  int ret = 1;

  if(strcmp(command, "digest") == 0) {
    fprintf(stderr, "-- Running in silico digest --\n");
    if(fasta_file == NULL) {
      fprintf(stderr, "FASTA file required (-f)\n");
      return 1;
    }
    if(restriction_seq == NULL) {
      fprintf(stderr, "Restriction sequence is required (-r)\n");
      return 1;
    } else {
      fprintf(stderr, "recognition seq is '%s'\n", restriction_seq);
      if(strcmp(restriction_seq, "DLE1") == 0 || strcmp(restriction_seq, "DLE-1") == 0) {
        fprintf(stderr, "setting recognition seq to CTTAAG\n");
        restriction_seq = "CTTAAG";
      }
    }
    // make a list of restriction seqs - that's what digest wants
    char** rseqs = malloc(1 * sizeof(char*));
    rseqs[0] = restriction_seq;
    c = digest_fasta(fasta_file, rseqs, 1);
    ret = write_cmap(&c, stdout);
  }

  if(strcmp(command, "label") == 0) {
    fprintf(stderr, "-- Running alignment-based labeling -> CMAP --\n");
    if(bam_file == NULL) {
      fprintf(stderr, "BAM file required (-a)\n");
      return 1;
    }
    c = get_cmap_from_bam(bam_file, covg_threshold);
    ret = write_cmap(&c, stdout);
  }

  else if(strcmp(command, "align") == 0 || strcmp(command, "dtw") == 0) {
    if(bnx_file == NULL) {
      fprintf(stderr, "BNX file (-b) required\n");
      return 1;
    }
    if(cmap_file == NULL) {
      fprintf(stderr, "CMAP file (-c) required\n");
      return 1;
    }
    FILE* o = stdout;

    fprintf(stderr, "# Loading '%s'...\n", bnx_file);
    time_t t0 = time(NULL);
    cmap b = read_bnx(bnx_file);
    time_t t1 = time(NULL);
    fprintf(stderr, "# Loaded %d molecules in %.2f seconds\n", b.n_maps, (t1-t0));
    if(end_mol == -1) {
      end_mol = b.n_maps;
    }

    fprintf(stderr, "# Loading '%s'...\n", cmap_file);
    t0 = time(NULL);
    cmap c = read_cmap(cmap_file);
    t1 = time(NULL);
    fprintf(stderr, "# Loaded CMAP '%s': %d maps w/%d recognition sites in %.2f seconds\n", cmap_file, c.n_maps, c.n_rec_seqs, t1-t0);

    int ret;
    if(strcmp(command, "align") == 0)
      ret = hash_cmap(b, c, o, q, chain_threshold, dtw_threshold, max_qgrams, read_limit, bin_size, min_frag, min_labels, start_mol, end_mol);
    else { // dtw

      int q, r, rv, a;
      result aln;
      for(q = start_mol; q <= end_mol; q++) {
        if(b.molecules[q].n_labels < min_labels) continue; // enforce minimum molecule labels
        result* alignments = malloc(c.n_maps * 2 * sizeof(result));
        uint32_t* qfrags = u32_get_fragments(b.molecules[q].labels, b.molecules[q].n_labels, 1, 0); // 1 is bin_size (no discretization)
        for(r = 0; r < c.n_maps; r++) {
          uint32_t* rfrags = u32_get_fragments(c.molecules[r].labels, c.molecules[r].n_labels, 1, 0);
          for(rv = 0; rv <= 1; rv++) {
            aln = dtw(qfrags, rfrags, b.molecules[q].n_labels, c.molecules[r].n_labels, -1, -1, 0.2, rv); // ins_score, del_score, neutral_deviation
            aln.ref = r;
            if(aln.failed) {
              fprintf(stderr, "Alignment failed of query %d to ref %d\n", q, r);
            }
            alignments[r + rv*c.n_maps] = aln;
          }
        }

        // sort alignments by (DTW) score decreasing
        ks_mergesort(aln_cmp, c.n_maps, alignments, 0);

        for(r = 0; r < c.n_maps*2; r++) {
          aln = alignments[r];
          if(aln.score < dtw_threshold) break;
          fprintf(o, "%u\t", b.molecules[q].id); // query id
          fprintf(o, "%u\t", c.molecules[aln.ref].id); // target id
          fprintf(o, "%u\t", aln.qrev); // query reverse?
          fprintf(o, "%u\t", aln.qstart); // query start idx
          fprintf(o, "%u\t", aln.qend); // query end idx
          fprintf(o, "%u\t", b.molecules[q].n_labels); // query len idx
          fprintf(o, "%u\t", b.molecules[q].labels[aln.qstart].position); // query start
          fprintf(o, "%u\t", b.molecules[q].labels[aln.qend > 0 ? aln.qend-1 : 0].position); // query end
          fprintf(o, "%u\t", b.molecules[q].length); // query len
          fprintf(o, "%u\t", aln.tstart); // ref start idx
          fprintf(o, "%u\t", aln.tend); // ref end idx
          fprintf(o, "%u\t", c.molecules[aln.ref].n_labels); // ref len idx
          fprintf(o, "%u\t", c.molecules[aln.ref].labels[aln.tstart].position); // ref start
          fprintf(o, "%u\t", c.molecules[aln.ref].labels[aln.tend > 0 ? aln.tend-1 : 0].position); // ref end
          fprintf(o, "%u\t", c.molecules[aln.ref].length); // ref len
          fprintf(o, "%f\t", aln.score); // dtw score
          // dtw path
          for(i = 0; i < kv_size(aln.path); i++) {
            fprintf(o, "%c", kv_A(aln.path, i) == 0 ? '.' : (kv_A(aln.path, i) == 1 ? 'I' : 'D'));
          }
          fprintf(o, "\n");
        }
        if(r == 0) {
          fprintf(o, "%u\t-\t-\t-\t-\t%u\t-\t-\t%u\t-\t-\t-\t-\t-\t-\t-\t-\n", b.molecules[q].id, b.molecules[q].n_labels, b.molecules[q].length);
        }
      }
      ret = 0;
    }

    // TODO: clean up cmap/bnx memory
  }

  else if(strcmp(command, "simulate") == 0) {
    if(fasta_file == NULL) {
      fprintf(stderr, "FASTA file required (-f)\n");
      return 1;
    }
    if(restriction_seq == NULL) {
      fprintf(stderr, "Restriction sequence is required (-r)\n");
      return 1;
    } else {
      if(strcmp(restriction_seq, "DLE1") == 0 || strcmp(restriction_seq, "DLE-1") == 0) {
        restriction_seq = "CTTAAG";
      }
    }
    if(coverage < FLT_EPSILON) {
      fprintf(stderr, "Coverage is required (-x)\n");
      return 1;
    }
    // make a list of restriction seqs - that's what digest wants
    char** rseqs = malloc(1 * sizeof(char*));
    rseqs[0] = restriction_seq;

    fprintf(stderr, "-- Running optical mapping simulation --\n");
    c = simulate_bnx(fasta_file, rseqs, 1, break_rate, fn, fp, stretch_mean, stretch_std, min_frag, coverage);
    fprintf(stderr, "Done simulating, writing to BNX...\n");
    ret = write_bnx(&c, stdout);

    if(source_outfile != NULL) {
      fprintf(stderr, "Writing truth/source positions to '%s'\n", source_outfile);
	    FILE *fp = fopen(source_outfile, "w");
      fprintf(fp, "ref_id\tstart_pos\n");
      for(i = 0; i < kv_size(c.source); i++) {
        fprintf(fp, "%u\t%u\n", kv_A(c.source, i).ref_id, kv_A(c.source, i).pos);
      }
      fclose(fp);
    }
  }

  // free everything
  return ret;
}
