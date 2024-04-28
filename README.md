# Puzzle Hi-C: an accurate scaffolding software

High-quality, chromosome-scale genomes are essential for genomic analyses. Analyses, including 3D genomics, epigenetics, and comparative genomics, rely on a high-quality genome assembly, which is often assembled with the assistance of Hi-C data. However, current Hi-C assisted assembling algorithms either generate ordering and orientation errors, or fail to assemble high-quality chromosome-level scaffolds. Here, we offer Puzzle Hi-C, which is software that uses Hi-C reads to assign accurately contigs or scaffolds to chromosomes. Puzzle Hi-C uses the triangle region instead of the square region to count interactions in a Hi-C heatmap. This strategy dramatically diminishes scaffolding interference caused by long-range interactions. This software also introduces a dynamic, triangle window strategy during assembling. The triangle window is initially small and expands with interactions to produce more effective clustering. We show that Puzzle Hi-C outperforms state-of-the-art tools for scaffolding.

## Installation
#### python 3.9.0
* biopython==1.81
* h5py==3.11.0
* networkx==3.2.1
* numpy==1.24.4
* pandas==2.2.2
* scipy==1.13.0
* matplotlib==3.8.2

#### Juicer 1.5.6
Installation of Juicer please refer to [Juicer](https://github.com/aidenlab/juicer)

### Install required python packages
```bash
pip install -r requirements.txt
```



## Usage
```bash
usage: main.py [-h] -c CLUSTERS [-p PREFIX] [-s BINSIZE] -m MATRIX -f FASTA
               [-t CUTOFF] [-i INIT_TRIANGLESIZE] [-n NCPUS] -j JUICER_TOOLS

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        Chromosomes number.
  -p PREFIX, --prefix PREFIX
                        Prefix
  -s BINSIZE, --binsize BINSIZE
                        The bin size. 
  -m MATRIX, --matrix MATRIX
                        The matrix file path.eg: merge_nodup.txt
  -f FASTA, --fasta FASTA
                        Scaffold fasta file.
  -t CUTOFF, --cutoff CUTOFF
                        Score cutoff.
  -i INIT_TRIANGLESIZE, --init_trianglesize INIT_TRIANGLESIZE
                        init_trianglesize.
  -n NCPUS, --ncpus NCPUS
                        Number of threads used for computering.
  -j JUICER_TOOLS, --juicer_tools JUICER_TOOLS
                        juicer_tools path.
                        
eg: python3 /public/home/lgl/software/main.py -c 5 -p Arabidopsis -s 10000 -t 0.35 -i 6 -m merged_nodups.txt -f ./ref/Arabidopsis_1M.fasta -j /public/home/lgl/software/juicer/PBS/scripts/juicer_tools -n 35

```

## Quick Start
1. Use Juicer to generate ```merged_nodups.txt```. Please refer to [Juicer](https://github.com/aidenlab/juicer).
2. Run puzzle Hi-C script.

## Example
```shell

```



