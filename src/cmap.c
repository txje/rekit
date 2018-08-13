#include <stdlib.h>
#include <stdio.h>

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

int write_cmap(cmap *c, char* fn) {
  FILE* fp = fopen(fn, "w");
  if(!fp) {
    fprintf(stderr, "Failed to write to file '%s'\n", fn);
    return 1;
  }

  size_t i, j, k;

  fprintf(fp, "# CMAP File Version:\t0.1\n");
  fprintf(fp, "# Label Channels:\t1\n");
  for(j = 0; j < c->n_rec_seqs; j++) {
    fprintf(fp, "# Nickase Recognition Site %u:\t%s\n", j+1, c->rec_seqs[j]);
  }
  fprintf(fp, "# Number of Consensus Nanomaps:\t%u\n", c->nmaps);
  fprintf(fp, "#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\tPosition\tStdDev\tCoverage\tOccurrence\n");
  fprintf(fp, "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint\n");

  for(i = 0; i < c->nmaps; i++) {
    for(j = 0; j < c->n_rec_seqs; j++) {
      for(k = 0; k < c->map_lengths[i][j]; k++) {
        fprintf(fp, "%u\t%f\t%u\t%u\t%u\t%f\t%f\t%u\t%u\n", i+1, c->ref_lengths[i], c->map_lengths[i][j], k, j+1, c->labels[i][j][k].position, c->labels[i][j][k].stdev, c->labels[i][j][k].coverage, c->labels[i][j][k].occurrence);
      }
    }
  }

  fclose(fp);
}
