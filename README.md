# Puzzle Hi-C: an accurate scaffolding software

High-quality, chromosome-scale genomes are essential for genomic analyses. Analyses, including 3D genomics, epigenetics, and comparative genomics, rely on a high-quality genome assembly, which is often assembled with the assistance of Hi-C data. However, current Hi-C assisted assembling algorithms either generate ordering and orientation errors, or fail to assemble high-quality chromosome-level scaffolds. Here, we offer Puzzle Hi-C, which is software that uses Hi-C reads to assign accurately contigs or scaffolds to chromosomes. Puzzle Hi-C uses the triangle region instead of the square region to count interactions in a Hi-C heatmap. This strategy dramatically diminishes scaffolding interference caused by long-range interactions. This software also introduces a dynamic, triangle window strategy during assembling. The triangle window is initially small and expands with interactions to produce more effective clustering. We show that Puzzle Hi-C outperforms state-of-the-art tools for scaffolding.

## Usage
```bash
usage: Puzzle_Hi-C_test_c [-h] -c CLUSTERS [-p PREFIX] [-s BINSIZE] -m MATRIX
                          -f FASTA [-t CUTOFF] [-i INIT_TRIANGLESIZE]
                          [-n NCPUS]

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        Chromosomes number!
  -p PREFIX, --prefix PREFIX
                        Prefix!
  -s BINSIZE, --binsize BINSIZE
                        The bin size!
  -m MATRIX, --matrix MATRIX
                        The matrix file path!
  -f FASTA, --fasta FASTA
                        Scaffold fasta file!
  -t CUTOFF, --cutoff CUTOFF
                        Score cutoff!
  -i INIT_TRIANGLESIZE, --init_trianglesize INIT_TRIANGLESIZE
                        init_trianglesize!
  -n NCPUS, --ncpus NCPUS
                        Number of threads used for computering!

```
## Install required packages
```bash
pip install -r requirements.txt
```
