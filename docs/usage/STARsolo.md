# STARsolo.smk

Snakemake workflow for gene count quantification from single-cell/nucleus RNA-seq data by [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/Sika-Zheng-Lab/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![STARsolo.smk rulegraph](https://github.com/Sika-Zheng-Lab/SnakeNgs/blob/develop/img/STARsolo_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Alignment using [STAR](https://github.com/alexdobin/STAR) with the parameter `--soloType CB_UMI_Simple`.
2. Make index for bam files using [samtools](http://www.htslib.org/).

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/STARsolo.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
star_index: /path/to/star_index
barcode_whitelist: /path/to/barcode_whitelist.txt
soloUMIlen: 12 # UMI length
soloCellFilter: EmptyDrops_CR # cell filtering type and parameters
clipAdapterType: CellRanger4 # adapter clipping type
outFilterScoreMin: 30 # minimum alignment score
soloCBmatchWLtype: 1MM_multi_Nbase_pseudocounts # matching the Cell Barcodes to the WhiteList
soloUMIfiltering: MultiGeneUMI_CR # type of UMI filtering
soloUMIdedup: 1MM_CR # type of UMI deduplication (collapsing) algorithm
soloBarcodeReadLength: 150 # length of the barcode read
```

- `experiment_table.tsv` should contain the following information:

``` text
sample  R1  R2
sample1 path/to/sample1_L001_R1.fastq.gz,path/to/sample1_L002_R1.fastq.gz path/to/sample1_L001_R2.fastq.gz,path/to/sample1_L002_R2.fastq.gz
sample2 path/to/sample2_L001_R1.fastq.gz,path/to/sample2_L002_R1.fastq.gz path/to/sample2_L001_R2.fastq.gz,path/to/sample2_L002_R2.fastq.gz
sample3 path/to/sample3_L001_R1.fastq.gz,path/to/sample3_L002_R1.fastq.gz path/to/sample3_L001_R2.fastq.gz,path/to/sample3_L002_R2.fastq.gz
```

`R1` and `R2` are comma-separated paths to the FASTQ files for read 1 and read 2, respectively.

- `path/to/star_index` is the directory containing the [STAR](https://github.com/alexdobin/STAR) index.

- `path/to/barcode_whitelist.txt` is the barcode whitelist file for the single-cell/nucleus RNA-seq data.

Please refer to the [STARsolo manual](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) for parameter details.

## Docker image used in the workflow

- [quay.io/biocontainers/star:2.7.11a--h0033a41_0](https://quay.io/repository/biocontainers/star)
- [quay.io/biocontainers/samtools:1.18--h50ea8bc_1](https://quay.io/repository/biocontainers/samtools)
