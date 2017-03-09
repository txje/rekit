rekit
-----
A toolkit for efficiently processing restriction maps and such


Tools may eventually include, in order:

  * .bnx and .cmap format parsers
  * in silico digestion -> rmap
  * in silico restriction mapping w/some error profile
  * pairwise rmap overlap
  * rmap overlap graph -> consensus map (cmap)
  * cmap -> reference alignment
  * structural variant prediction from cmap -> ref alignment

All code is distributed under the MIT license.

BNX parsing and associated data structures are originally from https://github.com/yanlinlin82/bntools
but have been significantly modified.


Installation
------------

    git clone https://github.com/txje/rekit
    cd rekit/src
    git clone https://github.com/attractivechaos/klib
    cd ..
    make


Usage
-----

    Usage: rekit [command] [options]
    Commands:
      ovl: compute MinHash/pairwise Jaccard similarity
      aln: compute dynamic time warping glocal (overlap) alignments
    Options:
      ovl <bnx> <q> <h> <seed> <threshold> <max_qgram_hits>
        bnx: A single BNX file containing rmaps
        q: Size of q-gram/k-mer to hash
        h: Number of hash functions to apply
        seed: Seed to random number generator
        threshold: Minimum number of k-mers to declare a match
        max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored
      aln <bnx> <threshold>
        bnx: A single BNX file containing rmaps
        threshold: Score threshold to report alignment
