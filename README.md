[![GitHub License](https://img.shields.io/github/license/NaotoKubota/SnakeNgs)](https://github.com/NaotoKubota/SnakeNgs/blob/main/LICENSE)
[![GitHub Release](https://img.shields.io/github/v/release/NaotoKubota/SnakeNgs?style=flat)](https://github.com/NaotoKubota/SnakeNgs/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/NaotoKubota/SnakeNgs)](https://github.com/NaotoKubota/SnakeNgs/releases)
[![Create Release and Build Docker Image](https://github.com/NaotoKubota/SnakeNgs/actions/workflows/release-docker-build-push.yaml/badge.svg)](https://github.com/NaotoKubota/SnakeNgs/actions/workflows/release-docker-build-push.yaml)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/snakengs)](https://hub.docker.com/r/naotokubota/snakengs)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/snakengs)](https://hub.docker.com/r/naotokubota/snakengs)

# SnakeNgs (v0.1.4)

Pipelines for NGS data analysis written in [Snakemake](https://snakemake.readthedocs.io/en/stable/). Each pipeline is designed to be executed with [Singularity](https://sylabs.io/singularity/) containers.

Please refer to the [documentation](https://naotokubota.github.io/SnakeNgs/) for more information.

## Contents

- Fetching
	- [ngsFetch](https://naotokubota.github.io/SnakeNgs/usage/ngsFetch/)
- RNA-seq
	- [preprocessing_RNAseq.smk](https://naotokubota.github.io/SnakeNgs/usage/preprocessing_RNAseq)
	- [preprocessing_RNAseq_single.smk](https://naotokubota.github.io/SnakeNgs/usage/preprocessing_RNAseq_single)
	- [rMATS.smk](https://naotokubota.github.io/SnakeNgs/usage/rMATS)
	- [SUPPA2_diffSplice.smk](https://naotokubota.github.io/SnakeNgs/usage/SUPPA2_diffSplice)
	- [Whippet.smk](https://naotokubota.github.io/SnakeNgs/usage/Whippet)
	- [MAJIQ.smk](https://naotokubota.github.io/SnakeNgs/usage/MAJIQ)
	- [LeafCutter.smk](https://naotokubota.github.io/SnakeNgs/usage/LeafCutter)
- snRNA-seq
	- [kb-nac.smk](https://naotokubota.github.io/SnakeNgs/usage/kb-nac)
- ChIP-seq
	- [preprocessing_ChIPseq.smk](https://naotokubota.github.io/SnakeNgs/usage/preprocessing_ChIPseq)
	- [preprocessing_ChIPseq_single.smk](https://naotokubota.github.io/SnakeNgs/usage/preprocessing_ChIPseq_single)
	- [callpeak_ChIPseq.smk](https://naotokubota.github.io/SnakeNgs/usage/callpeak_ChIPseq)
- ATAC-seq
	- [callpeak_ATACseq.smk](https://naotokubota.github.io/SnakeNgs/usage/callpeak_ATACseq)
	- [differential_ATACseq.smk](https://naotokubota.github.io/SnakeNgs/usage/differential_ATACseq)
	- [footprinting_ATACseq.smk](https://naotokubota.github.io/SnakeNgs/usage/footprinting_ATACseq)
- iCLIP-seq & HITSCLIP
	- [preprocessing_iCLIPseq.smk](https://naotokubota.github.io/SnakeNgs/usage/preprocessing_iCLIPseq)
	- [preprocessing_HITSCLIP.smk](https://naotokubota.github.io/SnakeNgs/usage/preprocessing_HITSCLIP)
- File format conversion
	- [bam2cram.smk](https://naotokubota.github.io/SnakeNgs/usage/bam2cram)