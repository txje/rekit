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
  if(c->n_maps <= 0) {
    fprintf(stderr, "Number of molecules header line not found or 0\n");
    return 1;
  }
	return 0;
}

int read_bnx_molecule(FILE *fp, cmap *c) {
  char** parts;
	char buf[16384];
	char token_buf[256];
  int i;
  const char* delim = "\t";
  
  // ------ line 0 ------
  if(fgets(buf, sizeof(buf), fp) == NULL) {
    return 1;
  }

  char* ln = strndup(buf, sizeof(buf));
  parts = malloc(sizeof(char*) * 20); // always exactly 20 fields in a 0 header line (for version 1.3)
  i = 0;
  char* token;
  token = strtok(ln, delim);
  while(token != NULL) {
    parts[i++] = token;
    token = strtok(NULL, delim);
  }
  int nfields = i;
  assert(nfields == 20);

  assert(strcmp(parts[0], "0") == 0);
  int mapid = atoi(parts[1]) - 1;
  assert(mapid < c->n_maps);

  c->ref_lengths[mapid] = (size_t)atof(parts[2]);
  c->map_lengths[mapid] = (size_t)atoi(parts[5]) + 1; // actual count is always 1 higher... maybe doesn't count end marker
  c->labels[mapid] = malloc(c->map_lengths[mapid] * sizeof(label));
  
  // ------ line 1 ------
  if(fgets(buf, sizeof(buf), fp) == NULL) {
    return 1;
  }
  token = strtok(buf, delim); // labelchannel
  int channel = atoi(token);
  token = strtok(NULL, delim); // first label pos
  i = 0;
  while(token != NULL) {
    assert(i < c->map_lengths[mapid]);
    c->labels[mapid][i].position = (uint32_t)atof(token);
    c->labels[mapid][i].channel = (uint8_t)channel;
    c->labels[mapid][i].occurrence = 0; // this is not used for molecule data, so it will always be 0
    c->labels[mapid][i].stdev = 0; // these will be set explicitly for all except the last label
    c->labels[mapid][i].coverage = 0;
    token = strtok(NULL, delim);
    i++;
  }
  // sscanf implicitly skips the newline

  // ------ line 2 ------
  if(fgets(buf, sizeof(buf), fp) == NULL) {
    return 1;
  }
  token = strtok(buf, delim); // labelchannel
  assert(strcmp(token, "QX11") == 0);
  token = strtok(NULL, delim); // first label pos
  i = 0;
  while(token != NULL) {
    assert(i < c->map_lengths[mapid]);
    c->labels[mapid][i].stdev = (uint32_t)atof(token); // THIS IS SNR, NOT stdev
    token = strtok(NULL, delim);
    i++;
  }

  // ------ line 3 ------
  if(fgets(buf, sizeof(buf), fp) == NULL) {
    return 1;
  }
  token = strtok(buf, delim); // labelchannel
  assert(strcmp(token, "QX12") == 0);
  token = strtok(NULL, delim); // first label pos
  i = 0;
  while(token != NULL) {
    assert(i < c->map_lengths[mapid]);
    c->labels[mapid][i].coverage = (uint32_t)atof(token); // THIS IS Intensity, NOT coverage
    token = strtok(NULL, delim);
    i++;
  }

  free(ln);
  free(parts);
  free(token);

	//next_line(fp, buf, sizeof(buf));
	return 0;
}

int write_bnx(cmap *c, FILE* fp) {
  //FILE* fp = fopen(fn, "w");
  if(!fp) {
    fprintf(stderr, "Failed to write cmap (invalid file pointer)\n");
    return 1;
  }

  size_t i, j, k;

  fprintf(fp, "# BNX File Version:\t1.3\n");
  fprintf(fp, "# Label Channels:\t%d\n", c->n_rec_seqs);
  for(j = 0; j < c->n_rec_seqs; j++) {
    fprintf(fp, "# Nickase Recognition Site %u:\t%s\n", j+1, c->rec_seqs[j]);
  }
  fprintf(fp, "# Bases per Pixel:\t%u\n", 500);
  fprintf(fp, "# Number of Molecules:\t%u\n", c->n_maps);
  fprintf(fp, "# Min Label SNR:\t%.2f\n", 0);
  fprintf(fp, "#0h LabelChannel  MoleculeID  Length  AvgIntensity  SNR NumberofLabels  OriginalMoleculeId  ScanNumber  ScanDirection ChipId  Flowcell  RunId Column  StartFOV  StartX  StartY  EndFOV  EndX  EndY  GlobalScanNumber\n");
  fprintf(fp, "#0f int  int   float  float float int int int int string  int int int int int int int int int int\n");
  fprintf(fp, "#1h LabelChannel  LabelPositions[N]\n");
  fprintf(fp, "#1f int float\n");
  fprintf(fp, "#Qh QualityScoreID  QualityScores[N]\n");
  fprintf(fp, "#Qf string  float[N]\n");
  fprintf(fp, "# Quality Score QX11: Label SNR for channel 1\n");
  fprintf(fp, "# Quality Score QX12: Label Intensity for channel 1\n");

  // ScanNumber is always 1, ScanDirection is unknown (-1), GlobalScanNumber is always 1, RunId is always 1
  for(i = 0; i < c->n_maps; i++) {
    fprintf(fp, "%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\tsim\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 0, i+1, (float)c->ref_lengths[i], 0.0, 0.0, c->map_lengths[i]-1, i+1, 1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1);
    fprintf(fp, "1");
    for(k = 0; k < c->map_lengths[i]; k++) {
      fprintf(fp, "\t%.2f", (float)c->labels[i][k].position);
    }
    // qualities have one fewer than lengths because the end position has a position but no quality
    fprintf(fp, "\nQX11");
    for(k = 0; k < c->map_lengths[i] - 1; k++) {
      fprintf(fp, "\t%.2f", (float)c->labels[i][k].stdev);
    }
    fprintf(fp, "\nQX12");
    for(k = 0; k < c->map_lengths[i] - 1; k++) {
      fprintf(fp, "\t%.2f", (float)c->labels[i][k].coverage);
    }
    fprintf(fp, "\n");
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
