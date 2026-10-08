# Puzzle Hi-C: an accurate scaffolding software

High-quality, chromosome-scale genomes are essential for genomic analyses. Analyses, including 3D genomics, epigenetics, and comparative genomics, rely on a high-quality genome assembly, which is often assembled with the assistance of Hi-C data. However, current Hi-C assisted assembling algorithms either generate ordering and orientation errors, or fail to assemble high-quality chromosome-level scaffolds. Here, we offer Puzzle Hi-C, which is software that uses Hi-C reads to assign accurately contigs or scaffolds to chromosomes. Puzzle Hi-C uses the triangle region instead of the square region to count interactions in a Hi-C heatmap. This strategy dramatically diminishes scaffolding interference caused by long-range interactions. This software also introduces a dynamic, triangle window strategy during assembling. The triangle window is initially small and expands with interactions to produce more effective clustering. We show that Puzzle Hi-C outperforms state-of-the-art tools for scaffolding.

## Recent updates
The latest update focuses on reducing memory usage and improving runtime efficiency while preserving the original scaffolding behavior.

Main improvements:

* Replaced the dense pair-orientation counting matrix with sparse accumulation for observed Hi-C links.
* Reduced the orientation matrix memory footprint by storing orientation values with a smaller integer type.
* Avoided large temporary upper-triangle index arrays during orientation matrix mirroring.
* Reused read-only worker context during multiprocessing on Linux to reduce repeated deserialization and duplicated worker memory.
* Kept compatibility fallbacks for environments that cannot inherit multiprocessing context.

In an _Arabidopsis thaliana_ 1 Mb contig benchmark on 35 CPU workers, peak memory usage was reduced from `19103504K` in the previous version to `6333308K` in the optimized version, while the test completed successfully and generated the expected `.agp`, `.fa`, `.Chrom.sizes`, and `.hic` outputs.

## Installation
#### python 3.9.0
* biopython==1.81
* h5py==3.11.0
* networkx==3.2.1
* numpy==1.24.4
* pandas==1.5.3
* scipy==1.13.0
* matplotlib==3.8.2

#### Optional alignment and visualization tools
* For BAM/SAM input: `pip install -r requirements-alignment.txt` (pysam).
* For `.hic` export: [Juicer tools](https://github.com/aidenlab/juicer). The Juicer
  alignment pipeline is not required for other input formats. Use `--skip-hic`
  to scaffold without Juicer tools or Java.

### Install required python packages
```bash
pip install -r requirements.txt  
```



## Usage
```bash
usage: main.py [-h] -c CLUSTERS -m MATRIX (-j JUICER_TOOLS | --skip-hic) -f FASTA [-p PREFIX] [-s BINSIZE] [-t CUTOFF] [-i INIT_TRIANGLESIZE]
               [-n NCPUS] [-e] [-g GAP]
               [--input-format {auto,juicer,bam,sam,pairs,validpairs}] [--min-mapq MIN_MAPQ]

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        Chromosomes number.
  -m MATRIX, --matrix MATRIX
                        Read-level contacts: Juicer, BAM/SAM, .pairs or HiC-Pro allValidPairs.
  -j JUICER_TOOLS, --juicer_tools JUICER_TOOLS
                        juicer_tools path for .hic export.
  --skip-hic            Skip .hic export; no Juicer installation is needed.
  --input-format {auto,juicer,bam,sam,pairs,validpairs}
                        Infer from filename by default; unknown extensions mean Juicer.
  --min-mapq MIN_MAPQ   Minimum MAPQ for both ends (default: 0).
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
1. Prepare deduplicated, filtered read-level Hi-C contacts against the same contigs
   as your FASTA, in any supported format below.
2. Run Puzzle Hi-C, choosing `--skip-hic` or `-j /path/to/juicer_tools`.
   Existing Juicer commands remain supported.

### Supported inputs (no Juicer alignment pipeline required)

`-m/--matrix` is retained for compatibility, but it takes **read-level contacts**,
not an aggregate matrix. Every accepted row/pair represents one contact.

* **Juicer** `merged_nodups.txt`: eight-column short or long read-pair format.
  Existing strand, 1-based position and restriction-fragment fields are retained.
  Nine-column weighted short format is rejected rather than silently losing counts.
* **BAM/SAM**: query-name sorted, paired alignments. Install the optional pysam
  dependency and name-sort first (an index is not needed):
  ```bash
  samtools sort -n -o reads.name.bam reads.bam
  python main.py -c 5 -m reads.name.bam -f contigs.fa --min-mapq 30 --skip-hic
  ```
  The header must declare `SO:queryname`; coordinate-sorted or merely collated
  files are rejected. Use unique template names, including across merged libraries.
  Both primary mates must be present. Orphans, ambiguous primary groups,
  unmapped, duplicate-marked and QC-failed pairs are excluded. Secondary records
  are ignored. Proper-pair status is **not** required, so inter-contig links are
  retained. Pairs with supplementary alignments or an `SA` tag are conservatively
  excluded: this adapter does not resolve chimeric Hi-C alignments. For those
  reads, use a Hi-C-aware parser such as [pairtools](https://pairtools.readthedocs.io/)
  and supply its filtered, deduplicated `.pairs` output instead.
* **4DN/pairtools `.pairs`**: standard seven columns
  `readID chr1 pos1 chr2 pos2 strand1 strand2`, or a `#columns:` header naming
  those columns in any order, with optional extra columns. `chrom1`/`chrom2`
  aliases are accepted. Without a header, exactly seven columns are required.
  Unmapped `!` ends are skipped. When `pair_type` is present, only `UU`, `UR` and
  `RU` contacts are retained; duplicate-marked `DD` and other types are skipped.
  ```bash
  python main.py -c 5 -m sample.pairs.gz -f contigs.fa --skip-hic
  ```
* **HiC-Pro `validPairs` / `allValidPairs`**: first seven columns
  `readID chr1 pos1 strand1 chr2 pos2 strand2`; later fragment and quality
  columns are allowed. Prefer filtered, deduplicated `allValidPairs` output:
  ```bash
  python main.py -c 5 -m sample.allValidPairs -f contigs.fa --skip-hic
  ```

Text inputs can be gzip-compressed (`.gz`). Blank lines and `#` comments are
ignored. Auto-detection recognizes `.bam`, `.sam`, `.pairs`, `.validPairs` and
`.allValidPairs` (case-insensitive, optionally followed by `.gz`); use
`--input-format pairs` or `--input-format validpairs` for other filenames.
Unknown extensions default to Juicer for backward compatibility. Malformed
records stop with a diagnostic instead of producing partial normalized input.

Positions are **1-based mapped 5′ endpoints**, never fabricated bin centers.
For BAM/SAM, the forward endpoint is the first aligned reference base; the reverse
endpoint is the last aligned reference base, including reference-consuming CIGAR
operations and excluding clipping. `.pairs` and HiC-Pro positions are retained.
Converted contacts use dummy Juicer fragment IDs 0 and 1 to avoid same-fragment
exclusion during `.hic` export; these are not inferred restriction sites.

`--min-mapq` defaults to 0 (no additional MAPQ filter). A positive threshold is
applied to **both** ends, and MAPQ 255 (unavailable) is excluded. Text input must
contain the quality fields: long Juicer columns 9/12, HiC-Pro columns 11/12, or
`.pairs` header columns `mapq1`/`mapq2`. Missing quality fields are an error when
filtering is requested. This adapter does not call duplicates or perform other
Hi-C protocol filtering; preprocess upstream and use the same reference FASTA.

**HiC-Pro `.matrix`, `.cool`, `.mcool` and `.hic` inputs are not supported.**
They aggregate contacts into bins and lose the individual endpoint coordinates
needed by Puzzle Hi-C's repeated scaffold rearrangement, triangle scoring and
error correction. Expanding counts at bin centers would change that algorithm;
normalized matrix weights are not molecule counts. Use HiC-Pro `allValidPairs`
instead. Supporting aggregate matrices requires a separate bin-aware algorithm.
CRAM is also not supported by this adapter.

`--skip-hic` keeps the normal `.agp`, `.fa` and `.Chrom.sizes` assembly outputs
but omits `.hic` generation. To export a `.hic`, replace it with
`-j /path/to/juicer_tools`; do not specify both options. Run each assembly in a
separate working directory, as the pipeline uses fixed intermediate filenames.

Format references: [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf),
[4DN pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md),
[pairtools](https://pairtools.readthedocs.io/en/latest/formats.html),
[HiC-Pro outputs](https://nservant.github.io/HiC-Pro/RESULTS.html),
[Juicer Pre](https://github.com/aidenlab/juicer/wiki/Pre#file-format).

### Tests

```bash
# Text adapters require only the Python standard library.
python -m unittest discover -s tests -v
# Install requirements.txt and requirements-alignment.txt to run BAM/SAM and
# full-pipeline integration tests too; otherwise optional tests report skips.
```

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
        -j ${Juicer}/PBS/scripts/juicer_tools \
        -n 8
```
## Utils
### Generate  ```.assembly``` and ```.hic``` files for Juicebox Assembly Tools (JBAT)
```shell
python agp2assembly.py -a agpfile -m merge_nodup.txt -j {Juicer}/PBS/scripts/juicer_tools -p prefix
```
### Generate  ```.agp``` file according to ```.assembly``` file
```shell
python agp2assembly.py -a assembly_file -g 100 -p prefix
```

## Citation
If you use Puzzle Hi-C in your work, please cite:

Lin G, Huang Z, Yue T, Chai J, Li Y, Yang H, et al. (2024) Puzzle Hi-C: An accurate scaffolding software. PLoS ONE 19(7): e0298564. https://doi.org/10.1371/journal.pone.0298564

## To do list

*  ~~Generate  ```.assembly``` file for Juicebox Assembly Tools (JBAT)~~
