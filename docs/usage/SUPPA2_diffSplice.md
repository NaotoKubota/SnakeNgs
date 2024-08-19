# SUPPA2_diffSplice.smk

Snakemake workflow for differential RNA splicing analysis using [SUPPA2](https://github.com/comprna/SUPPA).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![SUPPA2_diffSplice.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/SUPPA2_diffSplice_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Quantify transcript abundances (TPM) for each sample using [salmon](https://combine-lab.github.io/salmon/).
2. Merge the quantification results into a single file.
3. Make a list of alternative splicing events using [SUPPA2](https://github.com/comprna/SUPPA) `generateEvents`.
4. Quantify Percent Spliced In (PSI) values using [SUPPA2](https://github.com/comprna/SUPPA) `psiPerEvent`.
5. Separate the PSI values into two groups.
6. Separate the TPM values into two groups.
7. Differential splicing analysis using [SUPPA2](https://github.com/comprna/SUPPA) `diffSplice`.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/SUPPA2_diffSplice.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
salmon_index: /path/to/salmon_index
gtf: /path/to/reference_transcriptome.gtf
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

- `/path/to/salmon_index` is the directory containing the [salmon](https://combine-lab.github.io/salmon/) index.
- `/path/to/reference_transcriptome.gtf` is the reference transcriptome in GTF format (e.g. `Homo_sapiens.GRCh38.106.gtf` for human transcriptome).

## Docker image used in the workflow

- [combinelab/salmon:1.10.1](https://hub.docker.com/r/combinelab/salmon)
- [naotokubota/suppa:2.3](https://hub.docker.com/r/naotokubota/suppa)