# kb-nac.smk

Snakemake workflow for UMI counting from single-nucleus RNA-seq data with [kb-python](https://github.com/pachterlab/kb_python), a Python wrapper for [kallisto and bustools](https://www.kallistobus.tools/).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![kb-nac.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/kb-nac_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Build a kallisto index using the reference genome and transcriptome by `kb ref` with the parameter `--workflow nac`.
2. Quantify UMI counts using `kb count` with the parameter `--workflow nac --h5ad --gene-names --sum total --filter bustools --overwrite`.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/kb-nac.smk \
--configfile <path to config.yaml> \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
technology: 10xv3 # [10xv1, 10xv2, 10xv3, 10XV3_ULTIMA, BDWTA, INDROPSV3, Visium]
dna_fasta: /path/to/reference_dna.fa
gtf: /path/to/reference_transcriptome.gtf
```

- `experiment_table.tsv` should contain the following information:

``` text
sample  R1  R2
sample1 path/to/sample1_L001_R1.fastq.gz,path/to/sample1_L002_R1.fastq.gz path/to/sample1_L001_R2.fastq.gz,path/to/sample1_L002_R2.fastq.gz
sample2 path/to/sample2_L001_R1.fastq.gz,path/to/sample2_L002_R1.fastq.gz path/to/sample2_L001_R2.fastq.gz,path/to/sample2_L002_R2.fastq.gz
sample3 path/to/sample3_L001_R1.fastq.gz,path/to/sample3_L002_R1.fastq.gz path/to/sample3_L001_R2.fastq.gz,path/to/sample3_L002_R2.fastq.gz
```

`R1` and `R2` are comma-separated paths to the FASTQ files for read 1 and read 2, respectively.

- `technology` refer to the assay that generated the data. Supported assays are found in `kb --list`:

``` bash
$ singularity exec docker://quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2 kb --list
List of supported single-cell technologies

Positions syntax: `input file index, start position, end position`
When start & end positions are None, refers to the entire file
Custom technologies may be defined by providing a kallisto-supported technology string
(see https://pachterlab.github.io/kallisto/manual)

name            description                            on-list    barcode                    umi        cDNA
------------    -----------------------------------    -------    -----------------------    -------    -----------------------
10XV1           10x version 1                          yes        0,0,14                     1,0,10     2,None,None
10XV2           10x version 2                          yes        0,0,16                     0,16,26    1,None,None
10XV3           10x version 3                          yes        0,0,16                     0,16,28    1,None,None
10XV3_ULTIMA    10x version 3 sequenced with Ultima    yes        0,22,38                    0,38,50    0,62,None
BDWTA           BD Rhapsody                            yes        0,0,9 0,21,30 0,43,52      0,52,60    1,None,None
BULK            Bulk (single or paired)                                                                 0,None,None 1,None,None
CELSEQ          CEL-Seq                                           0,0,8                      0,8,12     1,None,None
CELSEQ2         CEL-SEQ version 2                                 0,6,12                     0,0,6      1,None,None
DROPSEQ         DropSeq                                           0,0,12                     0,12,20    1,None,None
INDROPSV1       inDrops version 1                                 0,0,11 0,30,38             0,42,48    1,None,None
INDROPSV2       inDrops version 2                                 1,0,11 1,30,38             1,42,48    0,None,None
INDROPSV3       inDrops version 3                      yes        0,0,8 1,0,8                1,8,14     2,None,None
SCRUBSEQ        SCRB-Seq                                          0,0,6                      0,6,16     1,None,None
SMARTSEQ2       Smart-seq2  (single or paired)                                                          0,None,None 1,None,None
SMARTSEQ3       Smart-seq3                                                                   0,11,19    0,11,None 1,None,None
SPLIT-SEQ       SPLiT-seq                                         1,10,18 1,48,56 1,78,86    1,0,10     0,None,None
STORMSEQ        STORM-seq                                                                    1,0,8      0,None,None 1,14,None
SURECELL        SureCell for ddSEQ                                0,0,6 0,21,27 0,42,48      0,51,59    1,None,None
Visium          10x Visium                             yes        0,0,16                     0,16,28    1,None,None
```

- `reference_dna.fa` is the reference genome in FASTA format (e.g., `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`).

- `reference_transcriptome.gtf` is the reference transcriptome in GTF format (e.g., `Homo_sapiens.GRCh38.106.gtf`).

Please refer to the [tutorial](../tutorial/snRNAseq_count.md) for more information.

## Docker image used in the workflow

- [quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2](https://quay.io/repository/biocontainers/kb-python)
