# preprocessing_RNAseq_single.smk

Snakemake workflow for preprocessing **single**-end bulk RNA-seq data.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![preprocessing_RNAseq_single.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/preprocessing_RNAseq_single_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Quality control using [fastp](https://github.com/OpenGene/fastp) with the default parameters.
2. Alignment using [STAR](https://github.com/alexdobin/STAR) with the parameter `--outFilterMultimapNmax 1`.
3. Convert the SAM file to BAM file and sort using [samtools](http://www.htslib.org/).
4. Make bigWig files using [deepTools](https://deeptools.readthedocs.io/en/develop/) `bamCoverage`.
5. Make summary statistics using [MultiQC](https://multiqc.info/).

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/preprocessing_RNAseq_single.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: path/to/output
samples: ["SRRXXXXXX", "SRRYYYYYY", "SRRZZZZZZ"]
star_index: path/to/star_index
```

- `path/to/output` should contain `fastq` directory with the following structure:

``` bash
output/
└── fastq
    ├── SRRXXXXXX.fastq.gz
    ├── SRRYYYYYY.fastq.gz
    └── SRRZZZZZZ.fastq.gz
```

- `path/to/star_index` is the directory containing the [STAR](https://github.com/alexdobin/STAR) index.

Please refer to the [tutorial](../tutorial/RNAseq_preprocessing.md) for more information.

## Docker image used in the workflow

- [quay.io/biocontainers/fastp:0.23.4--hadf994f_2](https://quay.io/repository/biocontainers/fastp)
- [quay.io/biocontainers/star:2.7.11a--h0033a41_0](https://quay.io/repository/biocontainers/star)
- [quay.io/biocontainers/samtools:1.18--h50ea8bc_1](https://quay.io/repository/biocontainers/samtools)
- [quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1](https://quay.io/repository/biocontainers/deeptools)
- [multiqc/multiqc:latest](https://hub.docker.com/r/multiqc/multiqc)
