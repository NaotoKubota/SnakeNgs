# SnakeNgs documentation

Pipelines for NGS data analysis written in [Snakemake](https://snakemake.readthedocs.io/en/stable/). Each pipeline is designed to be executed with [Singularity](https://sylabs.io/singularity/) containers.

## Contents

- [Installation](installation.md)
- Tutorial
    - [Fetching public NGS data](tutorial/Fetching.md)
    <!-- - [QC and mapping for RNA-seq](tutorial/RNAseq_preprocessing.md)
    - [UMI count from single-nucleus RNA-seq](tutorial/snRNAseq_count.md)
    - [QC and mapping for ChIP-seq/ATAC-seq](tutorial/ChIP-ATAC_preprocessing.md)
    - [Peak calling for ChIP-seq/ATAC-seq](tutorial/ChIP-ATAC_callpeak.md)
    - [Differential analysis for ATAC-seq](tutorial/ATACseq_diff.md)
    - [Footprinting for ATAC-seq](tutorial/ATACseq_footprinting.md)
- Usage
    - Fetching
        - [ngsFetch](usage/ngsFetch.md)
    - RNA-seq
        - [preprocessing_RNAseq.smk](usage/preprocessing_RNAseq.md)
        - [preprocessing_RNAseq_single.smk](usage/preprocessing_RNAseq_single.md)
        - [rMATS.smk](usage/rMATS.md)
        - [SUPPA2_diffSplice.smk](usage/SUPPA2_diffSplice.md)
        - [Whippet.smk](usage/Whippet.md)
        - [MAJIQ.smk](usage/MAJIQ.md)
        - [LeafCutter.smk](usage/LeafCutter.md)
    - snRNA-seq
        - [kb-nac.smk](usage/kb-nac.md)
    - ChIP-seq
        - [preprocessing_ChIPseq.smk](usage/preprocessing_ChIPseq.md)
        - [preprocessing_ChIPseq_single.smk](usage/preprocessing_ChIPseq_single.md)
        - [callpeak_ChIPseq.smk](usage/callpeak_ChIPseq.md)
    - ATAC-seq
        - [callpeak_ATACseq.smk](usage/callpeak_ATACseq.md)
        - [differential_ATACseq.smk](usage/differential_ATACseq.md)
        - [footprinting_ATACseq.smk](usage/footprinting_ATACseq.md)
    - iCLIP-seq & HITSCLIP
        - [preprocessing_iCLIPseq.smk](usage/preprocessing_iCLIPseq.md)
        - [preprocessing_HITSCLIP.smk](usage/preprocessing_HITSCLIP.md)
    - File format conversion
        - [bam2cram.smk](usage/bam2cram.md)
