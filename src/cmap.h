#include <stdint.h>

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

typedef struct label {
  uint32_t position;
  float stdev;
  uint16_t coverage;
  uint8_t channel;
  uint16_t occurrence;
} label;

typedef struct cmap {
  label** labels; // nested list by ref, then position
  size_t n_maps;
  size_t* ref_lengths;
  size_t* map_lengths;
  char** rec_seqs;
  size_t n_rec_seqs;
} cmap;

int write_cmap(cmap *c, FILE* fp);
cmap read_cmap(const char* fn);
