# Whippet.smk

Snakemake workflow for differential RNA splicing analysis using [Whippet](https://github.com/timbitz/Whippet.jl).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![Whippet.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/Whippet_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Quantify Percent Spliced In (PSI) values using [Whippet](https://github.com/timbitz/Whippet.jl) `quant`.
2. Differential splicing analysis using [Whippet](https://github.com/timbitz/Whippet.jl) `delta`.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/Whippet.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
whippet_index: /path/to/whippet_index/graph.jls
```

- `/path/to/output` is the directory where the output files will be saved.
- `/path/to/experiment_table.tsv` is a tab-separated file containing the following information:

``` text
sample	fastq	group
sample1	/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz	Ref
sample2	/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz	Ref
sample3	/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz	Ref
sample4	/path/to/sample4_R1.fastq.gz,/path/to/sample4_R2.fastq.gz	Alt
sample5	/path/to/sample5_R1.fastq.gz,/path/to/sample5_R2.fastq.gz	Alt
sample6	/path/to/sample6_R1.fastq.gz,/path/to/sample6_R2.fastq.gz	Alt
```

The `fastq` column should contain the paths to the FASTQ files for each sample. R1 and R2 files should be separated by a comma. The `group` column must be specified as `Ref` or `Alt`. This workflow will perform the differential splicing analysis between the two groups.

- `/path/to/whippet_index/graph.jls` is the [Whippet](https://github.com/timbitz/Whippet.jl) index file.

## Docker image used in the workflow

- [cloxd/whippet:1.6.1](https://hub.docker.com/r/cloxd/whippet)
