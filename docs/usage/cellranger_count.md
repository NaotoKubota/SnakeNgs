# cellranger_count.smk

Snakemake workflow for gene count quantification from single-cell/nucleus RNA-seq data by [CellRanger](https://github.com/10XGenomics/cellranger).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![cellranger_count.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/cellranger_count_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Make aliases for the input FASTQ files and convert the file names to the 10x Genomics format (e.g., `sample_S1_L001_R1_001.fastq.gz`).
2. Run `cellranger count` to quantify gene expression from the input FASTQ files.
3. Make summary statistics using [MultiQC](https://multiqc.info/).

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/cellranger_count.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
transcriptome: /path/to/reference_transcriptome
create_bam: "false" # ["true", "false"]
```

- `experiment_table.tsv` should contain the following information:

``` text
sample  R1  R2
sample1 path/to/sample1_L001_R1.fastq.gz,path/to/sample1_L002_R1.fastq.gz path/to/sample1_L001_R2.fastq.gz,path/to/sample1_L002_R2.fastq.gz
sample2 path/to/sample2_L001_R1.fastq.gz,path/to/sample2_L002_R1.fastq.gz path/to/sample2_L001_R2.fastq.gz,path/to/sample2_L002_R2.fastq.gz
sample3 path/to/sample3_L001_R1.fastq.gz,path/to/sample3_L002_R1.fastq.gz path/to/sample3_L001_R2.fastq.gz,path/to/sample3_L002_R2.fastq.gz
```

`R1` and `R2` are comma-separated paths to the FASTQ files for read 1 and read 2, respectively.

- `reference_transcriptome` should be the path to the reference transcriptome directory made by `cellranger mkref`.

- `create_bam` should be either `true` or `false`. If `true`, the workflow will create BAM files from the CellRanger output.

## Docker image used in the workflow

- [litd/docker-cellranger:v9.0.0](https://hub.docker.com/r/litd/docker-cellranger)
- [multiqc/multiqc:v1.27](https://hub.docker.com/r/multiqc/multiqc)
