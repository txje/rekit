rekit
-----
A toolkit for efficiently processing optical maps and such


Tools may eventually include, in order:

Existing features:
  * .bnx and .cmap format parsers
  * in silico digestion -> cmap
  * alignment-based labeling (BAM) -> cmap
  * in silico optical mapping (simulation) w/some error profile
  * feature-based molecule/cmap alignment with DTW refinement

Future features:
  * constructing consensus maps from pairwise molecule alignments
  * structural variant prediction from alignments

All code is distributed under the MIT license.

Some BNX parsing and associated data structures were borrowed from https://github.com/yanlinlin82/bntools
but have been significantly modified.

Requirements
------------

These need to be anywhere that the header and library files can be found by your compiler

  * zlib (https://zlib.net/)
  * htslib (https://github.com/samtools/htslib)

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
      overlap:  compute MinHash/pairwise Jaccard similarity
      align:    compute dynamic time warping glocal (overlap) alignments
      hash:     compute full q-gram intersection using a hash table
      simulate: simulate molecules
      digest:   in silico digestion
      label:    produce alignment-based reference CMAP
    Options:
      overlap  -b
      align    -b
      hash     -b
      simulate -frx --break-rate --fn --fp --min-frag --stretch-mean --stretch-std --source-output
      digest   -fr
      label    -a
        -b: bnx: A single BNX file containing molecules
        -c: cmap: A single CMAP file
        -f: fasta: Reference sequence to simulate from
        -a: bam: BAM alignment file
        -r: cutseq: Recognition/label site sequence
        -q: Size of q-gram/k-mer to hash (default: 4)
        -h: Number of hash functions to apply
        -t: threshold: Minimum number of q-grams to declare a match
        -m: max_qgram_hits: Maximum occurrences of a q-gram before it is considered repetitive and ignored
        -d: threshold: Score threshold to report alignment
        -x: Simulated molecule coverage
      simulate options:
        --break-rate: Probability of genome fragmentation per locus (default: 0.000005)
        --fn: Probability of missed label at true restriction site (default: 0.1)
        --fp: Probability of false-positive label (default: 0.05)
        --stretch-mean: Fragment stretch mean (default: 1.0)
        --stretch-std: Fragment stretch standard deviation (default: 0.05)
        --min-frag: Minimum detectable fragment size (default: 500)
        -s, --source-output: Output the reference positions of the simulated molecules to the given file
      label options:
        --coverage-threshold: Read coverage required (in ~300bp window) to call a label site (default: 10)

Output
------

Tab-delimited text to stdout with the following fields:

  1. Query (molecule) ID
  2. Reference map ID
  3. Reversed? (query)
  4. Query labels start index
  5. Query labels end index
  6. Query labels length
  7. Query label start position
  8. Query label end position
  9. Query total length
  10. Ref labels start index
  11. Ref labels end index
  12. Ref labels length
  13. Ref label start position
  14. Ref label end position
  15. Ref total length
  16. DTW alignment score
  17. DTW path string {'.', 'D', 'I'}
