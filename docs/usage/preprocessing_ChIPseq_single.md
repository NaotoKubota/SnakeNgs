# preprocessing_ChIPseq_single.smk

Snakemake workflow for preprocessing **single**-end ChIP-seq/ATAC-seq data.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![preprocessing_ChIPseq_single.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/preprocessing_ChIPseq_single_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Quality control using [fastp](https://github.com/OpenGene/fastp) with the default parameters.
2. Alignment using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) with the parameter specified in the `config.yaml`.
3. Convert the SAM file to BAM file and sort using [samtools](http://www.htslib.org/).
4. Remove duplicates using [Picard](https://broadinstitute.github.io/picard/) `MarkDuplicates` with the parameter `--REMOVE_DUPLICATES true`.
5. Make fingerprint plots using [deepTools](https://deeptools.readthedocs.io/en/develop/) `plotFingerprint` with the parameter `--minMappingQuality 30 --skipZeros --numberOfSamples 50000 --binSize 10000`.
6. Make bigWig files using [deepTools](https://deeptools.readthedocs.io/en/develop/) `bamCoverage` with the parameter `--binSize 1 --normalizeUsing CPM`.
7. Make summary statistics using [MultiQC](https://multiqc.info/).

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/preprocessing_ChIPseq_single.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: path/to/output
samples: ["SRRXXXXXX", "SRRYYYYYY", "SRRZZZZZZ"]
bowtie2_index: path/to/bowtie2_index
bowtie2_args: ""
```

- `path/to/output` should contain `fastq` directory with the following structure:

``` bash
output/
└── fastq
    ├── SRRXXXXXX.fastq.gz
    ├── SRRYYYYYY.fastq.gz
    └── SRRZZZZZZ.fastq.gz
```

- `path/to/bowtie2_index` is the directory containing the [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) index.

- `bowtie2_args` is the optional argument for Bowtie2 (e.g., `--very-sensitive`).

Please refer to the [tutorial](../tutorial/ChIP-ATAC_preprocessing.md) for more information.

## Docker image used in the workflow

- [quay.io/biocontainers/fastp:0.23.4--hadf994f_2](https://quay.io/repository/biocontainers/fastp)
- [quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0](https://quay.io/repository/biocontainers/bowtie2)
- [quay.io/biocontainers/samtools:1.18--h50ea8bc_1](https://quay.io/repository/biocontainers/samtools)
- [quay.io/biocontainers/picard:3.1.1--hdfd78af_0](https://quay.io/repository/biocontainers/picard)
- [quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1](https://quay.io/repository/biocontainers/deeptools)
- [multiqc/multiqc:v1.25](https://hub.docker.com/r/multiqc/multiqc)
