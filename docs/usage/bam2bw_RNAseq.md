# bam2bw_RNAseq.smk

Snakemake workflow for converting RNA-seq BAM files to BigWig files with [deepTools](https://deeptools.readthedocs.io/en/develop/) `bamCoverage`.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/Sika-Zheng-Lab/SnakeNgs) repository.

## Workflow

1. Read the experiment table (`sample`, `bam`) from the config file.
2. Convert each BAM file to BigWig using `bamCoverage` with `--binSize 1`.
3. Save all output BigWig files under `workdir/bigwig/`.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/bam2bw_RNAseq.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
```

`experiment_table.tsv` must include `sample` and `bam` columns:

``` text
sample	bam
sample1	/path/to/sample1.bam
sample2	/path/to/sample2.bam
```

- `/path/to/output` is the directory where logs, benchmarks, and BigWig files are saved.
- BigWig files are generated as `bigwig/{sample}.bw` in `workdir`.

## Docker image used in the workflow

- [quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1](https://quay.io/repository/biocontainers/deeptools)
