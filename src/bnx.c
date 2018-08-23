/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Linlin Yan, 2017 Jeremy Wang
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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include "cmap.h"

/*
...
# BNX File Version:  1.3
# Label Channels:  1
# Nickase Recognition Site 1:  cttaag;green_01
...
# Number of Molecules:  1322579
# Min Label SNR:  2.12
#0h LabelChannel  MoleculeID  Length  AvgIntensity  SNR NumberofLabels  OriginalMoleculeId  ScanNumber  ScanDirection ChipId  Flowcell  RunId Column  StartFOV  StartX  StartY  EndFOV  EndX  EndY  GlobalScanNumber
#0f int  int   float  float float int int int int string  int int int int int int int int int int
#1h LabelChannel  LabelPositions[N]
#1f int float
#Qh QualityScoreID  QualityScores[N]
#Qf string  float[N]
# Quality Score QX11: Label SNR for channel 1
# Quality Score QX12: Label Intensity for channel 1
0 1 937500.00 1315.06 128.58  147 1 1 -1  chips,SN_SN2ZXXGLPSMZFNWU,Run_08cc0975-68ad-4b04-ae4b-fd6ae1ef57a5,3553270258 1 1 5 2 341 1290  3 342 1963  1
1 2413.00 21025.00  24776.00  29189.00  32754.00  36917.00  39738.00  54746.00  55787.00  57873.00  67914.00  71401.00  73737.00  76457.00  87024.00...
QX11  30.34 56.27 24.36 33.02 34.97 46.89 86.29 39.59 58.28 44.01 35.01 39.02 67.79 32.06 ...
QX12  349.10  647.47  280.27  379.94  402.39  539.51  992.92  455.53  670.65  506.39  402.88 ...
*/

// it seems like "LabelChannel" is the first field and is always 0 for the header line and 1 (or higher, I've never seen it?) for the actual label line

int read_bnx_header(FILE *fp, cmap *c) {
	char buf[1024];
	int i, ch;
	char *val;
	while ((ch = getc(fp)) != EOF) {
		if (ch != '#') {
			ungetc(ch, fp);
			break;
		}
		if(fgets(buf, sizeof(buf), fp) == NULL) break; // returns a NULL pointer if there was nothing to read
		if(string_begins_with(buf, " BNX File Version:")) {
      assert(strcmp(get_val(buf), "1.3") == 0); // must be version 1.3
    }
		if(string_begins_with(buf, " Label Channels:")) {
      c->n_rec_seqs = atoi(get_val(buf));
      c->rec_seqs = malloc(c->n_rec_seqs * sizeof(char*));
		}
		if(string_begins_with(buf, " Nickase Recognition Site 1:")) {
      c->rec_seqs[0] = strdup(get_val(buf));
		}
		if(string_begins_with(buf, " Nickase Recognition Site 2:")) { // only supports up to 2 recog site right now
      c->rec_seqs[1] = strdup(get_val(buf));
		}
		if(string_begins_with(buf, " Number of Molecules:")) {
      c->n_maps = atoi(get_val(buf));
      c->map_lengths = malloc(c->n_maps * sizeof(size_t));
      c->ref_lengths = calloc(c->n_maps, sizeof(size_t));
      c->labels = malloc(c->n_maps * sizeof(label*));
    }
	}
	return 0;
}

int read_bnx_molecule(FILE *fp, cmap *c) {
  char** parts;
	char buf[1024];
  int i;
  const char* delim = "\t";
  
  // ------ line 0 ------
  if(fgets(buf, sizeof(buf), fp) == NULL) {
    return 1;
  }

  char* ln = strndup(buf, sizeof(buf));
  parts = malloc(sizeof(char*) * 21); // always exactly 21 fields in a 0 header line (for version 1.3)
  i = 0;
  char* token;
  token = strtok(ln, delim);
  while(token != NULL) {
    parts[i++] = token;
    token = strtok(NULL, delim);
  }
  int nfields = i;
  assert(nfields == 21);

  assert(strcmp(parts[0], "0") == 0);
  int mapid = atoi(parts[1]) - 1;
  assert(mapid < c->n_maps);

  c->ref_lengths[mapid] = (size_t)atof(parts[2]);
  c->map_lengths[mapid] = (size_t)atoi(parts[5]) + 1; // actual count is always 1 higher... maybe doesn't count end marker
  c->labels[mapid] = malloc(c->map_lengths[mapid] * sizeof(label));
  
  // ------ line 1 ------
  if(sscanf(buf, "%s", token) != 1) {
    fprintf(stderr, "Bad BNX format\n");
    return 1;
  }
  int channel = atoi(token);
  float f;
  for(i = 0; i < c->map_lengths[mapid]; i++) {
    if(sscanf(buf, "%f", &f) != 1) {
      fprintf(stderr, "Bad BNX format\n");
      return 1;
    }
    c->labels[mapid][i].position = (uint32_t)f;
    c->labels[mapid][site_id].channel = (uint8_t)channel;
    c->labels[mapid][site_id].occurrence = 0; // this is not used for molecule data, so it will always be 0
  }
  assert(fgetc(fp) == '\n');

  // ------ line 2 ------
  if(sscanf(buf, "%s", token) != 1) {
    fprintf(stderr, "Bad BNX format\n");
    return 1;
  }
  assert(strcmp(token, "QX11") == 0);
  float f;
  for(i = 0; i < c->map_lengths[mapid]; i++) {
    if(sscanf(buf, "%f", &f) != 1) {
      fprintf(stderr, "Bad BNX format\n");
      return 1;
    }
    c->labels[mapid][i].stdev = (uint32_t)f; // THIS IS SNR, NOT stdev
  }
  assert(fgetc(fp) == '\n');

  // ------ line 3 ------
  if(sscanf(buf, "%s", token) != 1) {
    fprintf(stderr, "Bad BNX format\n");
    return 1;
  }
  assert(strcmp(token, "QX12") == 0);
  float f;
  for(i = 0; i < c->map_lengths[mapid]; i++) {
    if(sscanf(buf, "%f", &f) != 1) {
      fprintf(stderr, "Bad BNX format\n");
      return 1;
    }
    c->labels[mapid][i].coverage = (uint32_t)f; // THIS IS Intensity, NOT coverage
  }
  assert(fgetc(fp) == '\n');

  free(ln);
  free(parts);
  free(token);

	//next_line(fp, buf, sizeof(buf));
	return 0;
}

int write_cmap(cmap *c, FILE* fp) {
  //FILE* fp = fopen(fn, "w");
  if(!fp) {
    fprintf(stderr, "Failed to write cmap (invalid file pointer)\n");
    return 1;
  }

  size_t i, j, k;

  fprintf(fp, "# CMAP File Version:\t0.1\n");
  fprintf(fp, "# Label Channels:\t%d\n", c->n_rec_seqs);
  for(j = 0; j < c->n_rec_seqs; j++) {
    fprintf(fp, "# Nickase Recognition Site %u:\t%s\n", j+1, c->rec_seqs[j]);
  }
  fprintf(fp, "# Number of Consensus Nanomaps:\t%u\n", c->n_maps);
  fprintf(fp, "#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n");
  fprintf(fp, "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n");

  for(i = 0; i < c->n_maps; i++) {
    for(k = 0; k < c->map_lengths[i]; k++) {
      // very weird - if the parameters are not cast, they can be arbitrarily reordered in the output string (presumably to optimize type-matching)
      fprintf(fp, "%u\t%.1f\t%u\t%u\t%u\t%.1f\t%.1f\t%u\t%u\n", i+1, (float)c->ref_lengths[i], c->map_lengths[i], k+1, c->labels[i][k].channel, (float)c->labels[i][k].position, (float)c->labels[i][k].stdev, c->labels[i][k].coverage, c->labels[i][k].occurrence);
    }
  }

  //fclose(fp);
  return 0;
}

cmap read_bnx(const char *filename) {
  cmap c;

	FILE *fp;

	fp = fopen(filename, "r");
	if (!fp) {
    fprintf(stderr, "File '%s' not found\n", filename);
		return c;
	}

	if (read_bnx_header(fp, &c) != 0) {
    fprintf(stderr, "File '%s' header could not be read\n", filename);
		fclose(fp);
		return c;
	}
	while (read_bnx_molecule(fp, &c) == 0) {
    // just looping through the file, reads 4 lines (1 molecule) each time
	}
	fclose(fp);

	return c;
}
