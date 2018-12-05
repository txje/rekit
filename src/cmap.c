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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include "cmap.h"

/*

# CMAP File Version:	0.1
# Label Channels:	1
# Nickase Recognition Site 1:	CTTAAG
# Number of Consensus Nanomaps:	66
#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
#f int	float	int	int	int	float	float	int	int
1	195471971.0	44559	1	1	3004108.0	1.0	1	1
1	195471971.0	44559	2	1	3010894.0	1.0	1	1
...

*/

int string_begins_with(char* s, char* pre) {
  return strlen(pre) <= strlen(s) && strncmp(pre, s, strlen(pre)) == 0;
}

// we stopped using this and are now assuming all lines will fit in the 1024-byte buffer BECAUSE we set the \n to \0 when we parse the value (below)
void next_line(FILE *fp, char *buf, size_t bufsize) {
	while (strchr(buf, '\n') == NULL) {
		if (fgets(buf, bufsize, fp) == NULL) break;
	}
}

char* get_val(char* buf) {
  char* val = strchr(buf, ':'); // pointer to ':'
  assert(val != NULL);
  ++val; // skip ':'
  while (*val && isblank(*val)) ++val; // skip whitespace
  val[strcspn(val, "\n")] = '\0'; // trim newline if present (otherwise just sets the null terminator to null again)
  return val;
}

int read_cmap_header(FILE *fp, cmap *c) {
	char buf[1024];
	int i, ch;
	char *val;
	while ((ch = getc(fp)) != EOF) {
		if (ch != '#') {
			ungetc(ch, fp);
			break;
		}
		if(fgets(buf, sizeof(buf), fp) == NULL) break; // returns a NULL pointer if there was nothing to read
    //fprintf(stderr, "%s\n", buf);
		if(string_begins_with(buf, " CMAP File Version:")) {
      assert(strcmp(get_val(buf), "0.1") == 0); // must be version 0.1
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
		if(string_begins_with(buf, " Number of Consensus Nanomaps:")) {
      c->n_maps = atoi(get_val(buf));
      c->molecules = malloc(c->n_maps * sizeof(molecule));
      for(i = 0; i < c->n_maps; i++) {
        c->molecules[i].length = 0;
      }
    }
		//next_line(fp, buf, sizeof(buf));
	}
	return 0;
}

static int read_cmap_line(FILE *fp, cmap *c) {
  char** parts;
	char buf[1024];
  int i;
  const char* delim = "\t";
  if(fgets(buf, sizeof(buf), fp) == NULL) {
    return 1;
  }

  char* ln = strndup(buf, sizeof(buf));
  parts = malloc(sizeof(char*) * 9); // always exactly 9 fields in a line (for version 0.1)
  i = 0;
  char* token;
  token = strtok(ln, delim);
  while(token != NULL) {
    parts[i++] = token;
    token = strtok(NULL, delim);
  }
  int nfields = i;
  assert(nfields == 9);

  int mapid = atoi(parts[0]) - 1;

  // if this is the first observed label for this ref/map, initialize it
  if(c->molecules[mapid].length == 0) {
    c->molecules[mapid].id = mapid+1;
    c->molecules[mapid].length = (size_t)atof(parts[1]);
    c->molecules[mapid].n_labels = (size_t)atoi(parts[2]) + 1; // actual count is always 1 higher... maybe doesn't count end marker
    c->molecules[mapid].labels = malloc(c->molecules[mapid].n_labels * sizeof(label));
  }
  int site_id = atoi(parts[3]) - 1;
  //printf("map %d, channel %d (%d labels), label %d\n", mapid, channel, c->map_lengths[mapid][channel], site_id);
  assert(site_id < c->molecules[mapid].n_labels);
  c->molecules[mapid].labels[site_id].position = (uint32_t)atof(parts[5]);
  c->molecules[mapid].labels[site_id].stdev = (float)atof(parts[6]);
  c->molecules[mapid].labels[site_id].coverage = (uint16_t)atoi(parts[7]);
  c->molecules[mapid].labels[site_id].channel = (uint8_t)atoi(parts[4]); // nickase channels start numbering at 1, 0 indicates end of the map/chromosome
  c->molecules[mapid].labels[site_id].occurrence = (uint16_t)atoi(parts[8]);

  free(ln);
  free(parts);

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
    for(k = 0; k < c->molecules[i].n_labels; k++) {
      // very weird - if the parameters are not cast, they can be arbitrarily reordered in the output string (presumably to optimize type-matching)
      fprintf(fp, "%u\t%.1f\t%u\t%u\t%u\t%.1f\t%.1f\t%u\t%u\n", c->molecules[i].id, (float)c->molecules[i].length, c->molecules[i].n_labels-1, k+1, c->molecules[i].labels[k].channel, (float)c->molecules[i].labels[k].position, (float)c->molecules[i].labels[k].stdev, c->molecules[i].labels[k].coverage, c->molecules[i].labels[k].occurrence);
    }
  }

  //fclose(fp);
  return 0;
}

cmap read_cmap(const char *fn) {
  cmap c;

  FILE *fp = fopen(fn, "r");
  if(!fp) {
    fprintf(stderr, "File '%s' not found\n", fn);
    return c;
  }

	if (read_cmap_header(fp, &c) != 0) {
    fprintf(stderr, "File '%s' header could not be read\n", fn);
		fclose(fp);
		return c;
	}
	while (read_cmap_line(fp, &c) == 0) {
    // we're just looping through all the lines
	}

	fclose(fp);
	return c;
}

void init_cmap(cmap* c) {
  c->n_maps = 0;
  c->molecules = NULL;
  c->rec_seqs = NULL;
  c->n_rec_seqs = 0;
  kv_init(c->source);
}

// positions should include the end pos of the chromosome
int add_map(cmap* c, uint32_t molid, uint32_t* positions, uint32_t n_pos, uint8_t channel) {
  uint32_t idx = c->n_maps;
  c->n_maps++;
  c->molecules = realloc(c->molecules, c->n_maps * sizeof(molecule));
  if(c->molecules == NULL) {
    fprintf(stderr, "Unable to allocate memory\n");
    return 1;
  };
  c->molecules[idx].id = molid;
  c->molecules[idx].length = positions[n_pos - 1];
  c->molecules[idx].n_labels = n_pos;
  c->molecules[idx].labels = malloc(n_pos * sizeof(label));
  uint32_t i;
  for(i = 0; i < n_pos; i++) {
    c->molecules[idx].labels[i].position = positions[i];
    if(i < n_pos - 1) {
      c->molecules[idx].labels[i].stdev = 1.0;
      c->molecules[idx].labels[i].coverage = 1;
      c->molecules[idx].labels[i].channel = channel;
      c->molecules[idx].labels[i].occurrence = 1;
    } else {
      c->molecules[idx].labels[i].stdev = 0.0;
      c->molecules[idx].labels[i].coverage = 1;
      c->molecules[idx].labels[i].channel = 0;
      c->molecules[idx].labels[i].occurrence = 0;
    }
  }
  return 0;
}

// filters labels based on minimum resolution (default: 500)
// filtered_labels should already have been alloc'd
size_t filter_labels(label* labels, size_t n_labels, label* filtered_labels, int resolution_min) {
  int i;
  size_t j = 0;
  for(i = 0; i < n_labels; i++) {
    if(1) { //j == 0 || labels[i].position - filtered_labels[j-1].position >= resolution_min) {
      filtered_labels[j++] = labels[i];
    } else {
      // sets the position of adjacent labels < min apart as the midpoint between them
      // DOES NOT account for other label features (incl. channel) - these are not used right now but they may be in the future
      filtered_labels[j-1].position = filtered_labels[j-1].position + (labels[i].position - filtered_labels[j-1].position) / 2;
    }
  }
  return j;
}
