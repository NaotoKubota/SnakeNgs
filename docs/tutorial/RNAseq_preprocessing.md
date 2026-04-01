# QC and mapping for RNA-seq

This tutorial walks you through the complete RNA-seq preprocessing workflow using the `preprocessing_RNAseq.smk` pipeline. You will download public paired-end RNA-seq data, prepare reference files, run the pipeline, and inspect the outputs.

By the end of this tutorial, you will have:

- Quality-controlled FASTQ files
- STAR-aligned BAM files
- BigWig coverage tracks
- RNA-seq QC metrics
- A MultiQC summary report

!!! note "Prerequisites"

    Before starting, make sure you have the following installed and configured:

    - [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) (≥ 3.7)
    - [Snakemake](https://snakemake.readthedocs.io/en/stable/) (≥ 7.0)
    - [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository cloned locally
    - [ngsfetch](https://github.com/NaotoKubota/ngsfetch) for downloading FASTQ files
    - ~20 GB of free disk space (for reference genome, STAR index, and output files)

## 1. Download example data

In this tutorial, we use three paired-end RNA-seq samples from a study on nonsense-mediated mRNA decay (NMD) in mouse cortical development ([GSE295221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295221), BioProject [PRJNA1253720](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1253720); Lin et al., *Cell Rep* 2026). These are total RNA-seq libraries (TruSeq Stranded Total RNA) from E18.5 mouse cortex.

| Accession | Sample | Genotype | Instrument | Layout | Read length |
|-----------|--------|----------|------------|--------|-------------|
| [SRR33238275](https://www.ncbi.nlm.nih.gov/sra/SRR33238275) | GSM8943899 (Control Rep1) | Upf2 fl/fl | NovaSeq X Plus | Paired-end | 2 × 150 bp |
| [SRR33238274](https://www.ncbi.nlm.nih.gov/sra/SRR33238274) | GSM8943900 (Control Rep2) | Upf2 fl/+ | NovaSeq X Plus | Paired-end | 2 × 150 bp |
| [SRR33238273](https://www.ncbi.nlm.nih.gov/sra/SRR33238273) | GSM8943901 (Control Rep3) | Upf2 fl/fl | NovaSeq X Plus | Paired-end | 2 × 150 bp |

Create a working directory and download the FASTQ files using [ngsfetch](https://github.com/NaotoKubota/ngsfetch):

``` bash
# Create working directory
mkdir -p rnaseq_tutorial
cd rnaseq_tutorial

# Download paired-end RNA-seq data from SRA
ngsfetch -i SRR33238275 -o fastq -p 16
ngsfetch -i SRR33238274 -o fastq -p 16
ngsfetch -i SRR33238273 -o fastq -p 16
```

After downloading, your directory should look like this:

``` bash
rnaseq_tutorial/
└── fastq
    ├── SRR33238273_1.fastq.gz
    ├── SRR33238273_2.fastq.gz
    ├── SRR33238274_1.fastq.gz
    ├── SRR33238274_2.fastq.gz
    ├── SRR33238275_1.fastq.gz
    └── SRR33238275_2.fastq.gz
```

## 2. Prepare reference files

The pipeline requires a STAR genome index and a GTF annotation file. Here we use the mouse reference genome (GRCm38/mm10) from Ensembl.

### Download the GTF file

``` bash
wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz
```

### Build the STAR index

!!! note

    Building a STAR index for the mouse genome requires ~32 GB of RAM and takes approximately 30–60 minutes. If you already have a STAR index, you can skip this step.

``` bash
# Download the genome FASTA
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Build STAR index
mkdir -p star_index
STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \
    --sjdbGTFfile Mus_musculus.GRCm38.102.gtf \
    --runThreadN 8
```

## 3. Create the configuration file

Create a `config.yaml` file in the working directory:

``` yaml
workdir: /path/to/rnaseq_tutorial
samples: ["SRR33238275", "SRR33238274", "SRR33238273"]
star_index: /path/to/star_index
gtf: /path/to/Mus_musculus.GRCm38.102.gtf
layout: "paired"
```

Replace `/path/to/` with the actual absolute paths on your system.

## 4. Run the pipeline

Execute the pipeline with Snakemake:

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/preprocessing_RNAseq.smk \
    --configfile config.yaml \
    --cores 8 \
    --use-singularity \
    --rerun-incomplete
```

!!! note

    Adjust `--cores` based on the number of CPU cores available on your system. The first run will pull the required Singularity containers automatically.

## 5. Inspect the outputs

Once the pipeline completes, the output directory will have the following structure:

``` bash
rnaseq_tutorial/
├── fastp/
│   ├── SRR33238275_1.fastq.gz
│   ├── SRR33238275_2.fastq.gz
│   ├── ...
│   └── log/
│       ├── SRR33238275.html
│       ├── SRR33238275.json
│       └── ...
├── star/
│   ├── SRR33238275/
│   │   ├── SRR33238275_Aligned.out.bam
│   │   └── SRR33238275_Log.final.out
│   └── ...
├── metrics/
│   ├── SRR33238275.picard.analysis.CollectRnaSeqMetrics
│   ├── SRR33238275.picard.analysis.CollectInsertSizeMetrics
│   └── ...
├── bigwig/
│   ├── SRR33238275.bw
│   └── ...
└── multiqc/
    └── multiqc_report.html
```

### Key output files

- **`fastp/log/*.html`** — Per-sample QC reports showing read quality, adapter content, and filtering statistics.
- **`star/*_Aligned.out.bam`** — Sorted BAM files with uniquely mapped reads (`--outFilterMultimapNmax 1`).
- **`star/*_Log.final.out`** — STAR alignment summary with mapping rates.
- **`metrics/*.CollectRnaSeqMetrics`** — Picard RNA-seq metrics including the percentage of reads mapping to coding, UTR, intronic, and intergenic regions.
- **`metrics/*.CollectInsertSizeMetrics`** — Insert size distribution for paired-end data.
- **`bigwig/*.bw`** — Normalized coverage tracks that can be loaded in genome browsers such as [IGV](https://igv.org/) or the [UCSC Genome Browser](https://genome.ucsc.edu/).
- **`multiqc/multiqc_report.html`** — Aggregated QC report combining fastp, STAR, and Picard metrics across all samples.

Open `multiqc/multiqc_report.html` in a web browser to review the overall quality of your experiment at a glance.

## 6. Summary and next steps

In this tutorial, you ran the `preprocessing_RNAseq.smk` pipeline to perform quality control, alignment, and metric collection for paired-end RNA-seq data.

For detailed parameter descriptions and single-end mode, see the [usage documentation](../usage/preprocessing_RNAseq.md).

The aligned BAM files produced by this pipeline can be used as input for downstream analyses available in SnakeNgs:

- [rMATS](../usage/rMATS.md)
- [SUPPA2](../usage/SUPPA2_diffSplice.md)
- [Whippet](../usage/Whippet.md)
- [MAJIQ](../usage/MAJIQ.md)
- [LeafCutter](../usage/LeafCutter.md)

The BAM files can also be used as input for [Shiba](https://github.com/Sika-Zheng-Lab/Shiba), a tool for differential RNA splicing analysis.
