# rMATS.smk

Snakemake workflow for differential RNA splicing analysis using [rMATS](http://rnaseq-mats.sourceforge.net/).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![rMATS.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/rMATS_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Parse the experiment table to create the sample list.
2. Run [rMATS](http://rnaseq-mats.sourceforge.net/) with `--novelSS` and the parameters specified in the `config.yaml` file.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/rMATS.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
gtf: /path/to/reference_transcriptome.gtf # Use unzipped GTF file

# rMATS
layout: paired # single or paired
readLength: 101 # Read length of the RNA-seq data
cstat: 0.1 # Cutoff for the statistical test
anchorLength: 6 # Anchor length for the reads
mil: 0 # Minimum intron length
mel: 10000 # Maximum exon length
```

- `/path/to/output` is the directory where the output files will be saved.
- `/path/to/experiment_table.tsv` is a tab-separated file, which is same as the one used in [Shiba](https://github.com/NaotoKubota/Shiba).

``` text
sample	bam	group
sample1	/path/to/sample1.bam	Ref
sample2	/path/to/sample2.bam	Ref
sample3	/path/to/sample3.bam	Ref
sample4	/path/to/sample4.bam	Alt
sample5	/path/to/sample5.bam	Alt
sample6	/path/to/sample6.bam	Alt
```

The `group` column must be specified as `Ref` or `Alt`. This workflow will perform the differential splicing analysis between the two groups.

- `/path/to/reference_transcriptome.gtf` is the reference transcriptome in GTF format (e.g. `Homo_sapiens.GRCh38.106.gtf` for human transcriptome). Note that the GTF file should be unzipped.

## Docker image used in the workflow

- [xinglab/rmats:v4.3.0](https://hub.docker.com/r/xinglab/rmats)
