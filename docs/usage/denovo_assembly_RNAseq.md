# denovo_assembly_RNAseq.smk

Snakemake workflow for **de novo transcriptome assembly** from short-read RNA-seq data using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq), with ORF prediction by [TransDecoder](https://github.com/TransDecoder/TransDecoder) and gene annotation by [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) against UniProt/Swiss-Prot. Supports both **paired-end** and **single-end** reads via the `layout` parameter.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

### Paired-end (`layout: "paired"`)

<figure markdown="span">
	![denovo_assembly_RNAseq.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/denovo_assembly_RNAseq_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

### Single-end (`layout: "single"`)

<figure markdown="span">
	![denovo_assembly_RNAseq_single.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/denovo_assembly_RNAseq_single_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

### Pipeline steps

1. Quality control using [fastp](https://github.com/OpenGene/fastp) with the default parameters.
2. De novo transcriptome assembly using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) with all samples combined.
3. Assembly statistics using Trinity's `TrinityStats.pl`.
4. ORF prediction using [TransDecoder](https://github.com/TransDecoder/TransDecoder) `LongOrfs`.
5. Homology search of assembled transcripts against UniProt/Swiss-Prot using [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) `blastx`.
6. Homology search of predicted peptides against UniProt/Swiss-Prot using BLAST+ `blastp`.
7. Coding region prediction using TransDecoder `Predict` with BLAST evidence.
8. Annotation summary table combining BLAST results and TransDecoder predictions.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/denovo_assembly_RNAseq.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--singularity-args "--bind $HOME:$HOME --bind /path/to/swissprot_db_dir" \
--rerun-incomplete
```

!!! warning

    Make sure to add the directory containing the Swiss-Prot BLAST database to `--singularity-args "--bind ..."` so that it is accessible from within the Singularity container.

`config.yaml` should contain the following information:

### Paired-end

``` yaml
workdir: path/to/output
samples: ["SRRXXXXXX", "SRRYYYYYY", "SRRZZZZZZ"]
layout: "paired"
max_memory: "32G"
swissprot_db: path/to/uniprot_sprot
```

- `path/to/output` should contain `fastq` directory with the following structure:

``` bash
output/
└── fastq
    ├── SRRXXXXXX_1.fastq.gz
    ├── SRRXXXXXX_2.fastq.gz
    ├── SRRYYYYYY_1.fastq.gz
    ├── SRRYYYYYY_2.fastq.gz
    ├── SRRZZZZZZ_1.fastq.gz
    └── SRRZZZZZZ_2.fastq.gz
```

### Single-end

``` yaml
workdir: path/to/output
samples: ["SRRXXXXXX", "SRRYYYYYY", "SRRZZZZZZ"]
layout: "single"
max_memory: "32G"
swissprot_db: path/to/uniprot_sprot
```

- `path/to/output` should contain `fastq` directory with the following structure:

``` bash
output/
└── fastq
    ├── SRRXXXXXX.fastq.gz
    ├── SRRYYYYYY.fastq.gz
    └── SRRZZZZZZ.fastq.gz
```

### Common settings

- `max_memory`: Maximum memory for Trinity (e.g. `"32G"`, `"64G"`). Trinity requires substantial memory; **32 GB or more is recommended**.
- `swissprot_db`: Path to the UniProt/Swiss-Prot BLAST database (without file extension). You need to prepare the database beforehand:

``` bash
# Download and build Swiss-Prot BLAST database
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot
```

## Output

| File | Description |
|------|-------------|
| `trinity/Trinity.fasta` | Assembled transcriptome |
| `trinity/Trinity_stats.txt` | Assembly statistics (N50, number of transcripts, etc.) |
| `transdecoder/Trinity.fasta.transdecoder.pep` | Predicted peptide sequences |
| `transdecoder/Trinity.fasta.transdecoder.bed` | Predicted coding regions in BED format |
| `transdecoder/Trinity.fasta.transdecoder.gff3` | Predicted coding regions in GFF3 format |
| `blast/blastx_swissprot.outfmt6` | BLASTx results (transcripts vs Swiss-Prot) |
| `blast/blastp_swissprot.outfmt6` | BLASTp results (predicted peptides vs Swiss-Prot) |
| `annotation/annotation_summary.tsv` | Annotation summary table |

## Docker images used in the workflow

- [quay.io/biocontainers/fastp:0.23.4--hadf994f_2](https://quay.io/repository/biocontainers/fastp)
- [trinityrnaseq/trinityrnaseq:2.15.2](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq)
- [quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0](https://quay.io/repository/biocontainers/transdecoder)
- [quay.io/biocontainers/blast:2.16.0--hc155240_2](https://quay.io/repository/biocontainers/blast)
