# QC, mapping, peak calling for ChIP-seq/ATAC-seq

This tutorial walks you through the complete ChIP-seq analysis workflow using two SnakeNgs pipelines: `preprocessing_ChIPseq.smk` for read preprocessing and `callpeak_ChIPseq.smk` for peak calling. You will download public ChIP-seq data, prepare reference files, run both pipelines sequentially, and inspect the outputs.

By the end of this tutorial, you will have:

- Quality-controlled and deduplicated BAM files
- Fingerprint plots for enrichment assessment
- BigWig coverage tracks
- MACS2 peak calls
- MultiQC summary reports

!!! note "Prerequisites"

    Before starting, make sure you have the following installed and configured:

    - [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) (≥ 3.7)
    - [Snakemake](https://snakemake.readthedocs.io/en/stable/) (≥ 7.0)
    - [SnakeNgs](https://github.com/Sika-Zheng-Lab/SnakeNgs) repository cloned locally
    - [ngsfetch](https://github.com/Sika-Zheng-Lab/ngsfetch) for downloading FASTQ files
    - ~20 GB of free disk space (for reference genome, Bowtie2 index, and output files)

## 1. Download example data

In this tutorial, we use three human ChIP-seq samples from a K562 cell line study ([PRJEB72776](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB72776)): two target ChIP samples (H3K4me3 and H3K27ac) sharing one input control. We download paired-end FASTQ files using [ngsfetch](https://github.com/Sika-Zheng-Lab/ngsfetch).

| Accession | Sample | Type | Instrument | Layout | Read length | Download size |
|-----------|--------|------|------------|--------|-------------|---------------|
| [ERR12721584](https://www.ncbi.nlm.nih.gov/sra/ERR12721584) | H3K4me3 | ChIP target | HiSeq 1500 | Paired-end | 2 × 51 bp | ~1.3 Gb |
| [ERR12721583](https://www.ncbi.nlm.nih.gov/sra/ERR12721583) | H3K27ac | ChIP target | HiSeq 1500 | Paired-end | 2 × 51 bp | ~1.3 Gb |
| [ERR12721582](https://www.ncbi.nlm.nih.gov/sra/ERR12721582) | Input (shared) | Input control | HiSeq 1500 | Paired-end | 2 × 51 bp | ~1.4 Gb |

``` bash
# Create working directory
mkdir -p chipseq_tutorial
cd chipseq_tutorial

# Download ChIP-seq data (2 targets + 1 input control)
ngsfetch -i ERR12721584 -o fastq -p 16
ngsfetch -i ERR12721583 -o fastq -p 16
ngsfetch -i ERR12721582 -o fastq -p 16
```

After downloading, your directory should look like this:

``` bash
chipseq_tutorial/
└── fastq
    ├── ERR12721582_1.fastq.gz
    ├── ERR12721582_2.fastq.gz
    ├── ERR12721583_1.fastq.gz
    ├── ERR12721583_2.fastq.gz
    ├── ERR12721584_1.fastq.gz
    └── ERR12721584_2.fastq.gz
```

In this example, `ERR12721584` (H3K4me3) and `ERR12721583` (H3K27ac) are target (IP) samples, and `ERR12721582` is the input control shared between both marks.

## 2. Prepare reference files

The preprocessing pipeline requires a Bowtie2 genome index for alignment.

### Download and build the Bowtie2 index

!!! note

    Building a Bowtie2 index for the human genome requires ~8 GB of RAM and takes approximately 30–60 minutes. If you already have a Bowtie2 index, you can skip this step.

``` bash
# Download the genome FASTA
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build Bowtie2 index
mkdir -p bowtie2_index
bowtie2-build --threads 8 \
    Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    bowtie2_index/GRCh38
```

## Step A: Preprocessing

The first step uses `preprocessing_ChIPseq.smk` to perform quality control, alignment, duplicate removal, and coverage track generation.

### 3. Create the preprocessing configuration file

Create a `config_preprocessing.yaml` file:

``` yaml
workdir: /path/to/chipseq_tutorial
samples: ["ERR12721584", "ERR12721583", "ERR12721582"]
bowtie2_index: /path/to/bowtie2_index/GRCh38
bowtie2_args: "--very-sensitive"
layout: "paired"
```

Replace `/path/to/` with the actual absolute paths on your system.

!!! note

    Include **all** samples (both target and input control) in the `samples` list for preprocessing, as all of them need to be aligned and deduplicated.

### 4. Run the preprocessing pipeline

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/preprocessing_ChIPseq.smk \
    --configfile config_preprocessing.yaml \
    --cores 8 \
    --use-singularity \
    --rerun-incomplete
```

### 5. Inspect the preprocessing outputs

Once the pipeline completes, the output directory will have the following structure:

``` bash
chipseq_tutorial/
├── fastp/
│   ├── ERR12721584_1.fq.gz
│   ├── ERR12721584_2.fq.gz
│   ├── ...
│   └── log/
│       ├── ERR12721584.html
│       └── ...
├── bowtie2/
│   ├── ERR12721584.sort.rmdup.bam
│   ├── ERR12721584.sort.rmdup.bam.bai
│   ├── ERR12721583.sort.rmdup.bam
│   ├── ERR12721583.sort.rmdup.bam.bai
│   ├── ERR12721582.sort.rmdup.bam
│   └── ERR12721582.sort.rmdup.bam.bai
├── bigwig/
│   ├── ERR12721584.bw
│   ├── ERR12721583.bw
│   └── ERR12721582.bw
├── plotFingerprint/
│   ├── fingerprint.png
│   ├── fingerprint.qc.txt
│   └── fingerprint.tab
├── metrics/
│   ├── ERR12721584.picard.analysis.CollectInsertSizeMetrics
│   └── ...
└── multiqc_preprocessing/
    └── multiqc_report.html
```

#### Key output files

- **`bowtie2/*.sort.rmdup.bam`** — Sorted and deduplicated BAM files ready for peak calling.
- **`plotFingerprint/fingerprint.png`** — Fingerprint plot generated by deepTools. This plot helps assess the enrichment quality of your ChIP-seq experiment—target samples should show a clear deviation from the diagonal compared to input.
- **`bigwig/*.bw`** — CPM-normalized coverage tracks for visualization in genome browsers.
- **`multiqc_preprocessing/multiqc_report.html`** — Aggregated QC report combining fastp, Bowtie2, and Picard metrics.

## Step B: Peak calling

The second step uses `callpeak_ChIPseq.smk` to call peaks with MACS2, comparing target samples against input control.

### 6. Create the experiment table

Create an `experiment_table.tsv` file that maps each target sample to its input control. Use **tab-separated** columns:

``` text
sample	target	control
ERR12721584	/path/to/chipseq_tutorial/bowtie2/ERR12721584.sort.rmdup.bam	/path/to/chipseq_tutorial/bowtie2/ERR12721582.sort.rmdup.bam
ERR12721583	/path/to/chipseq_tutorial/bowtie2/ERR12721583.sort.rmdup.bam	/path/to/chipseq_tutorial/bowtie2/ERR12721582.sort.rmdup.bam
```

Replace `/path/to/` with the actual absolute paths on your system. The `target` column points to the IP BAM file and the `control` column points to the input BAM file.

### 7. Create the peak calling configuration file

Create a `config_callpeak.yaml` file:

``` yaml
general:
  workdir: /path/to/chipseq_tutorial
  experiment_table: /path/to/chipseq_tutorial/experiment_table.tsv

macs2:
  broad: False
  broad_cutoff:
  genomesize: hs
  fileformat: BAMPE
  qvalue: 0.05

bedgraphtobigwig:
  assembly: hg38
```

Key parameters:

- `genomesize` — Effective genome size. Use `hs` for human, `mm` for mouse, `ce` for *C. elegans*, `dm` for *Drosophila*.
- `fileformat` — Use `BAMPE` for paired-end BAM files, or `BAM` for single-end.
- `broad` — Set to `True` for broad histone marks (e.g., H3K27me3, H3K36me3). Use `False` for narrow peaks (e.g., transcription factors, H3K4me3).

### 8. Run the peak calling pipeline

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/callpeak_ChIPseq.smk \
    --configfile config_callpeak.yaml \
    --cores 8 \
    --use-singularity \
    --rerun-incomplete
```

### 9. Inspect the peak calling outputs

``` bash
chipseq_tutorial/
├── macs2/
│   ├── ERR12721584/
│   │   ├── ERR12721584_peaks.xls
│   │   ├── ERR12721584_peaks.narrowPeak
│   │   ├── ERR12721584_summits.bed
│   │   ├── ERR12721584_treat_pileup.bdg
│   │   └── ERR12721584_treat_pileup.bw
│   └── ERR12721583/
│       ├── ERR12721583_peaks.xls
│       ├── ERR12721583_peaks.narrowPeak
│       ├── ERR12721583_summits.bed
│       ├── ERR12721583_treat_pileup.bdg
│       └── ERR12721583_treat_pileup.bw
└── multiqc_callpeak/
    └── multiqc_report.html
```

#### Key peak calling output files

- **`macs2/*_peaks.narrowPeak`** — Called peaks in ENCODE narrowPeak format with coordinates, significance, and fold enrichment.
- **`macs2/*_peaks.xls`** — Detailed peak information including p-values and q-values.
- **`macs2/*_summits.bed`** — Peak summit positions (single base-pair resolution).
- **`macs2/*_treat_pileup.bw`** — Treatment pileup converted to BigWig format for genome browser visualization.
- **`multiqc_callpeak/multiqc_report.html`** — Aggregated report for peak calling statistics.

## Note on ATAC-seq

The same preprocessing pipeline (`preprocessing_ChIPseq.smk`) can be used for ATAC-seq data. For peak calling on ATAC-seq, use the `callpeak_ATACseq.smk` pipeline instead of `callpeak_ChIPseq.smk`. The key differences are:

- ATAC-seq peak calling does **not** require a control (input) sample.
- The config uses a sample list instead of an experiment table.
- Use `fileformat: BAMPE` for paired-end ATAC-seq data.

For details, see the [callpeak_ATACseq.smk usage documentation](../usage/callpeak_ATACseq.md).

## 10. Summary and next steps

In this tutorial, you ran two pipelines sequentially:

1. `preprocessing_ChIPseq.smk` — QC, Bowtie2 alignment, deduplication, fingerprint plotting, and BigWig generation.
2. `callpeak_ChIPseq.smk` — MACS2 peak calling with input control normalization.

For detailed parameter descriptions, see the usage documentation for [preprocessing_ChIPseq.smk](../usage/preprocessing_ChIPseq.md) and [callpeak_ChIPseq.smk](../usage/callpeak_ChIPseq.md).

The peak files produced by these pipelines can be used as input for downstream analyses available in SnakeNgs:

- [differential_ATACseq.smk](../usage/differential_ATACseq.md) — Differential accessibility analysis with DESeq2
- [footprinting_ATACseq.smk](../usage/footprinting_ATACseq.md) — Transcription factor footprinting with TOBIAS
- [footprinting_timeseries_ATACseq.smk](../usage/footprinting_timeseries_ATACseq.md) — Time-series footprinting analysis
