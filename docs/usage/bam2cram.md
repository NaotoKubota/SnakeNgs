# bam2cram.smk

Snakemake workflow for file format conversion from BAM to CRAM to compress the file size.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![bam2cram.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/bam2cram_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Convert BAM files to CRAM files using [samtools](http://www.htslib.org/).
2. Make CRAM index files (`.crai`) using [samtools](http://www.htslib.org/).

## Usage

``` bash
snakemake -s /rhome/naotok/SnakeNgs/Snakefile_bam2cram \
--config \
workdir=/path/to/output \
reference=/path/to/reference_genome.fa \
bam_list=/path/to/bam_list.txt \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

- `path/to/output` is the directory where the log files will be saved.
- `path/to/reference_genome.fa` is the reference genome file in FASTA format (e.g., `Mus_musculus.GRCm38.dna_sm.primary_assembly.fa` for mouse genome).
- `path/to/bam_list.txt` is a text file containing the list of BAM files to be converted to CRAM files:

``` text
/path/to/SRRXXXXXX.bam
/path/to/SRRYYYYYY.bam
/path/to/SRRZZZZZZ.bam
/path/to/SRRAAAAAA.bam
/path/to/SRRBBBBBB.bam
/path/to/SRRCCCCCC.bam
```

The CRAM files will be saved in the same directory as the input BAM files with the `.cram` extension.

## Docker image used in the workflow

- [quay.io/biocontainers/samtools:1.18--h50ea8bc_1](https://quay.io/repository/biocontainers/samtools)
