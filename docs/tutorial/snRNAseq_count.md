# UMI count from single-nucleus RNA-seq

This tutorial walks you through the single-nucleus RNA-seq (snRNA-seq) UMI counting workflow using the `kb-nac.smk` pipeline. You will download a public 10x Chromium snRNA-seq dataset, prepare reference files, run the pipeline, and inspect the outputs.

By the end of this tutorial, you will have:

- A kallisto index
- Filtered and unfiltered count matrices in h5ad format
- BUS file inspection reports
- A MultiQC summary report

!!! note "Prerequisites"

    Before starting, make sure you have the following installed and configured:

    - [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) (≥ 3.7)
    - [Snakemake](https://snakemake.readthedocs.io/en/stable/) (≥ 7.0)
    - [SnakeNgs](https://github.com/Sika-Zheng-Lab/SnakeNgs) repository cloned locally
    - [ngsfetch](https://github.com/Sika-Zheng-Lab/ngsfetch) for downloading FASTQ files
    - ~20 GB of free disk space (for reference genome and output files)

## 1. Download example data

In this tutorial, we use two samples from a study on nonsense-mediated mRNA decay (NMD) in mouse cortical development ([GSE295222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295222), BioProject [PRJNA1253721](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1253721); Lin et al., *Cell Rep* 2026). This dataset profiled E17.5 mouse cortex nuclei using the 10x Genomics Chromium platform.

| Sample | Accessions (2 lanes) | Genotype | Instrument | Layout |
|--------|---------------------|----------|------------|--------|
| GSM8943907 (Control) | [SRR33238112](https://www.ncbi.nlm.nih.gov/sra/SRR33238112), [SRR33238113](https://www.ncbi.nlm.nih.gov/sra/SRR33238113) | Upf2 fl/+ | NovaSeq X Plus | Paired-end |
| GSM8943908 (Upf2cKO) | [SRR33238110](https://www.ncbi.nlm.nih.gov/sra/SRR33238110), [SRR33238111](https://www.ncbi.nlm.nih.gov/sra/SRR33238111) | Upf2 fl/fl;Emx1-Cre | NovaSeq X Plus | Paired-end |

Each sample was sequenced across two lanes on an Illumina NovaSeq X Plus.

Download the FASTQ files using [ngsfetch](https://github.com/Sika-Zheng-Lab/ngsfetch):

``` bash
# Create working directory
mkdir -p snrnaseq_tutorial
cd snrnaseq_tutorial

# Download snRNA-seq data from SRA (2 samples × 2 lanes)
ngsfetch -i SRR33238112 -o fastq -p 16
ngsfetch -i SRR33238113 -o fastq -p 16
ngsfetch -i SRR33238110 -o fastq -p 16
ngsfetch -i SRR33238111 -o fastq -p 16
```

!!! note

    For 10x Chromium data, `_1.fastq.gz` typically contains the cell barcode + UMI (R1) and `_2.fastq.gz` contains the cDNA insert (R2). Verify that your files follow this convention.

## 2. Prepare the experiment table

Create an `experiment_table.tsv` file that maps sample names to their FASTQ file paths:

``` text
sample	R1	R2
control	/path/to/snrnaseq_tutorial/fastq/SRR33238112_1.fastq.gz,/path/to/snrnaseq_tutorial/fastq/SRR33238113_1.fastq.gz	/path/to/snrnaseq_tutorial/fastq/SRR33238112_2.fastq.gz,/path/to/snrnaseq_tutorial/fastq/SRR33238113_2.fastq.gz
Upf2cKO	/path/to/snrnaseq_tutorial/fastq/SRR33238110_1.fastq.gz,/path/to/snrnaseq_tutorial/fastq/SRR33238111_1.fastq.gz	/path/to/snrnaseq_tutorial/fastq/SRR33238110_2.fastq.gz,/path/to/snrnaseq_tutorial/fastq/SRR33238111_2.fastq.gz
```

Replace `/path/to/` with the actual absolute paths on your system.

!!! note

    If your samples were sequenced across multiple lanes, provide comma-separated paths for each lane:

    ``` text
    sample1	path/to/sample1_L001_R1.fastq.gz,path/to/sample1_L002_R1.fastq.gz	path/to/sample1_L001_R2.fastq.gz,path/to/sample1_L002_R2.fastq.gz
    ```

## 3. Prepare reference files

The pipeline requires a reference genome FASTA and a GTF annotation file. The kallisto index will be built automatically by the pipeline using `kb ref`.

Download the mouse reference files from Ensembl:

``` bash
# Download the genome FASTA
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Download the GTF file
wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz
```

## 4. Create the configuration file

Create a `config.yaml` file in the working directory:

``` yaml
workdir: /path/to/snrnaseq_tutorial
experiment_table: /path/to/snrnaseq_tutorial/experiment_table.tsv
technology: 10xv3
dna_fasta: /path/to/Mus_musculus.GRCm38.dna.primary_assembly.fa
gtf: /path/to/Mus_musculus.GRCm38.102.gtf
```

Replace `/path/to/` with the actual absolute paths on your system.

The `technology` parameter specifies the single-cell assay. Common options include:

| Technology | Description |
|------------|-------------|
| `10xv2` | 10x Chromium v2 |
| `10xv3` | 10x Chromium v3 |
| `BDWTA` | BD Rhapsody |
| `INDROPSV3` | inDrops v3 |
| `Visium` | 10x Visium spatial |

For the full list of supported technologies, see the [kb-nac.smk usage documentation](../usage/kb-nac.md).

## 5. Run the pipeline

Execute the pipeline with Snakemake:

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/kb-nac.smk \
    --configfile config.yaml \
    --cores 8 \
    --use-singularity \
    --rerun-incomplete
```

!!! note

    The first step (`kb ref`) builds the kallisto index, which can take 30–60 minutes and requires ~16 GB of RAM for the mouse genome. The index is built once and reused for all samples.

## 6. Inspect the outputs

Once the pipeline completes, the output directory will have the following structure:

``` bash
snrnaseq_tutorial/
├── kb_index/
│   ├── index.idx
│   ├── t2g.txt
│   ├── cdna.fa
│   ├── intron.fa
│   ├── cdna_t2c.txt
│   └── intron_t2c.txt
├── kb/
│   ├── control/
│   │   ├── counts_unfiltered/
│   │   │   └── adata.h5ad
│   │   ├── counts_filtered/
│   │   │   └── adata.h5ad
│   │   └── inspect.json
│   └── Upf2cKO/
│       ├── counts_unfiltered/
│       │   └── adata.h5ad
│       ├── counts_filtered/
│       │   └── adata.h5ad
│       └── inspect.json
└── multiqc/
    └── multiqc_report.html
```

### Key output files

- **`kb_index/`** — Kallisto index files built by `kb ref`. These include the index, transcript-to-gene mapping, and cDNA/intron FASTA files.
- **`kb/*/counts_unfiltered/adata.h5ad`** — Unfiltered count matrix containing all barcodes, in [AnnData](https://anndata.readthedocs.io/) h5ad format.
- **`kb/*/counts_filtered/adata.h5ad`** — Filtered count matrix containing only cell-associated barcodes (filtered by bustools).
- **`kb/*/inspect.json`** — BUS file inspection summary with statistics on the number of reads, barcodes, and UMIs.
- **`multiqc/multiqc_report.html`** — Aggregated QC summary report.

### Loading the count matrix in Python

The h5ad output can be directly loaded with [Scanpy](https://scanpy.readthedocs.io/):

``` python
import scanpy as sc

adata = sc.read_h5ad("kb/control/counts_filtered/adata.h5ad")
print(adata)
```

## 7. Alternative tools

SnakeNgs also provides alternative pipelines for single-cell/nucleus RNA-seq quantification:

- [STARsolo](../usage/STARsolo.md) — Gene count quantification using STAR's built-in single-cell mode. Produces Cell Ranger-compatible output.
- [Cell Ranger](../usage/cellranger_count.md) — 10x Genomics' official pipeline for gene expression quantification.

## 8. Summary and next steps

In this tutorial, you ran the `kb-nac.smk` pipeline to build a kallisto index and quantify UMI counts from single-nucleus RNA-seq data.

For detailed parameter descriptions, see the [usage documentation](../usage/kb-nac.md).

The filtered count matrices can be used for downstream analysis with tools such as:

- [Scanpy](https://scanpy.readthedocs.io/) — Python toolkit for single-cell analysis
- [Seurat](https://satijalab.org/seurat/) — R toolkit for single-cell analysis
