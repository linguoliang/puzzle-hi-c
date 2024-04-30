# Puzzle Hi-C: an accurate scaffolding software

High-quality, chromosome-scale genomes are essential for genomic analyses. Analyses, including 3D genomics, epigenetics, and comparative genomics, rely on a high-quality genome assembly, which is often assembled with the assistance of Hi-C data. However, current Hi-C assisted assembling algorithms either generate ordering and orientation errors, or fail to assemble high-quality chromosome-level scaffolds. Here, we offer Puzzle Hi-C, which is software that uses Hi-C reads to assign accurately contigs or scaffolds to chromosomes. Puzzle Hi-C uses the triangle region instead of the square region to count interactions in a Hi-C heatmap. This strategy dramatically diminishes scaffolding interference caused by long-range interactions. This software also introduces a dynamic, triangle window strategy during assembling. The triangle window is initially small and expands with interactions to produce more effective clustering. We show that Puzzle Hi-C outperforms state-of-the-art tools for scaffolding.

## Installation
#### python 3.9.0
* biopython==1.81
* h5py==3.11.0
* networkx==3.2.1
* numpy==1.24.4
* pandas==1.5.3
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
usage: main.py [-h] -c CLUSTERS -m MATRIX -j JUICER_TOOLS -f FASTA [-p PREFIX] [-s BINSIZE] [-t CUTOFF] [-i INIT_TRIANGLESIZE]
               [-n NCPUS] [-e] [-g GAP]

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        Chromosomes number.
  -m MATRIX, --matrix MATRIX
                        The matrix file path.eg: merge_nodup.txt
  -j JUICER_TOOLS, --juicer_tools JUICER_TOOLS
                        juicer_tools path.
  -f FASTA, --fasta FASTA
                        Scaffold fasta file.
  -p PREFIX, --prefix PREFIX
                        Output prefix! Default: sample.
  -s BINSIZE, --binsize BINSIZE
                        The bin size. Default: 10000.
  -t CUTOFF, --cutoff CUTOFF
                        Score cutoff, 0.25-0.5 recommended. default: 0.3.
  -i INIT_TRIANGLESIZE, --init_trianglesize INIT_TRIANGLESIZE
                        Initial triangle size. Default: 3.
  -n NCPUS, --ncpus NCPUS
                        Number of threads. Default: 1.
  -e, --error_correction
                        For error correction! Default: False.
  -g GAP, --gap GAP     The size of gap between scaffolds. Default: 100.

                        
eg: python3 /public/home/lgl/software/puzzle-hic/main.py -c 5 -p Arabidopsis -s 10000 -t 0.35 -i 6 -m merged_nodups.txt -f ./ref/Arabidopsis_1M.fasta -j /public/home/lgl/software/juicer/PBS/scripts/juicer_tools -n 35

```

## Quick Start
1. Use Juicer to generate ```merged_nodups.txt```. Please refer to [Juicer](https://github.com/aidenlab/juicer).
2. Run puzzle Hi-C script.

## Example
We use _Arabidopsis thaliana_ as an example. The T2T genome we used is from  NCBI: [GCA_028009825.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028009825.2/). Put the genome in ```ref``` directory. We split this genome into 1Mb length contigs using [generate_test_data.py](utils%2Fgenerate_test_data.py).
The Hi-C data is downloaded from GSA: [CRR302669](https://ngdc.cncb.ac.cn/gsa/browse/CRA004538/CRR302669). Put the Hi-C data in ```raw_data``` directory. Now, here is a step-by-step guide.

```shell
# Set up juicer env
export Juicer=/path/to/juicer  # here is your Juicer path.
# Set up Puzzle Hi-C env
export PuzzleHiC=/path/to/Puzzle Hi-C  # here is your Puzzle Hi-C path.
# Split Genomes
cd ref
python3 ${PuzzleHiC}/generate_test_data.py GCA_028009825.2_Col-CC_genomic.fna Arabidopsis   # Output: Arabidopsis_1M.fasta
# Build bwa index
bwa index Arabidopsis_1M.fasta
# Create  restriction site file. 
python2 ${Juicer}/PBS/scripts/generate_site_positions.py DpnII Arabidopsis_1M Arabidopsis_1M.fasta
# Prepare for running Juicedr
cd ..
mkdir -p juicer_1M/fastq
cd juicer_1M/fastq
ln -s ../../raw_data/CRR302669_f1.fastq.gz CRR302669_R1.fastq.gz
ln -s ../../raw_data/CRR302669_r2.fastq.gz CRR302669_R2.fastq.gz

# Run Juicer
cd ..
${Juicer}/CPU/juicer.sh \
        -t 8 \
        -y ../ref/Arabidopsis_1M_DpnII.txt \
        -p ../ref/Arabidopsis_1M.chrom.sizes \
        -z ../ref/Arabidopsis_1M.fasta
# Run Puzzle Hi-C
cd ..
mkdir Puzzle_hic_1M
cd Puzzle_hic_1M
ln -s ../juicer_1M/aligned/merged_nodups.txt ./
python3 ${PuzzleHiC}/main.py \
        -c 5 -p Arabidopsis -s 10000 \
        -t 0.35 -i 6 -m merged_nodups.txt \
        -f ../ref/Arabidopsis_1M.fasta \
        -j {Juicer}/PBS/scripts/juicer_tools \
        -n 8
```



