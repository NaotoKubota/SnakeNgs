#!/usr/bin/env bash
# Generate dummy test data for Snakefile dry-run testing in CI
set -euo pipefail

BASE="$(cd "${GITHUB_WORKSPACE:-.}" && pwd)/tests/testdata"

rm -rf "$BASE"

###############################################################################
# Shared resources
###############################################################################

# STAR index (chrLength.txt needed at parse-time by functions.smk)
mkdir -p "$BASE/star_index"
printf '200000000\n180000000\n' > "$BASE/star_index/chrLength.txt"
printf 'chr1\t200000000\nchr2\t180000000\n' > "$BASE/star_index/chrNameLength.txt"
touch "$BASE/star_index/Genome" "$BASE/star_index/SA" "$BASE/star_index/SAindex"
touch "$BASE/star_index/genomeParameters.txt"

# Bowtie2 index (prefix: $BASE/bowtie2_index/genome)
mkdir -p "$BASE/bowtie2_index"
for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
    touch "$BASE/bowtie2_index/genome.$ext"
done

# Dummy annotation/reference files
touch "$BASE/dummy.gtf" "$BASE/dummy.gff3"
touch "$BASE/dummy.fa" "$BASE/dummy.fa.fai"
touch "$BASE/motifs.jaspar" "$BASE/barcode_whitelist.txt" "$BASE/whippet_index.jls"
mkdir -p "$BASE/salmon_index" "$BASE/transcriptome"

###############################################################################
# Dummy BAM / peak / fastq files (referenced by experiment tables)
###############################################################################

mkdir -p "$BASE/bams" "$BASE/peaks" "$BASE/fastqs" "$BASE/chipseq_bams"
for s in sample1 sample2 sample3 sample4 sample5 sample6; do
    touch "$BASE/bams/$s.bam" "$BASE/peaks/${s}_peaks.narrowPeak"
    touch "$BASE/fastqs/${s}_1.fastq.gz" "$BASE/fastqs/${s}_2.fastq.gz"
done
touch "$BASE/chipseq_bams/target1.bam" "$BASE/chipseq_bams/control1.bam"

###############################################################################
# Experiment tables
###############################################################################

mkdir -p "$BASE/experiment_tables"

# bam + peak + group (differential_ATACseq, footprinting_ATACseq)
{
    printf 'sample\tbam\tpeak\tgroup\n'
    for i in 1 2 3; do printf "sample%s\t%s/bams/sample%s.bam\t%s/peaks/sample%s_peaks.narrowPeak\tRef\n" "$i" "$BASE" "$i" "$BASE" "$i"; done
    for i in 4 5 6; do printf "sample%s\t%s/bams/sample%s.bam\t%s/peaks/sample%s_peaks.narrowPeak\tAlt\n" "$i" "$BASE" "$i" "$BASE" "$i"; done
} > "$BASE/experiment_tables/bam_peak_group.tsv"

# bam + peak + group timeseries (footprinting_timeseries_ATACseq)
{
    printf 'sample\tbam\tpeak\tgroup\n'
    printf "sample1\t%s/bams/sample1.bam\t%s/peaks/sample1_peaks.narrowPeak\tESC\n" "$BASE" "$BASE"
    printf "sample2\t%s/bams/sample2.bam\t%s/peaks/sample2_peaks.narrowPeak\tNPC\n" "$BASE" "$BASE"
    printf "sample3\t%s/bams/sample3.bam\t%s/peaks/sample3_peaks.narrowPeak\td25\n" "$BASE" "$BASE"
    printf "sample4\t%s/bams/sample4.bam\t%s/peaks/sample4_peaks.narrowPeak\td25\n" "$BASE" "$BASE"
    printf "sample5\t%s/bams/sample5.bam\t%s/peaks/sample5_peaks.narrowPeak\td50\n" "$BASE" "$BASE"
    printf "sample6\t%s/bams/sample6.bam\t%s/peaks/sample6_peaks.narrowPeak\td50\n" "$BASE" "$BASE"
} > "$BASE/experiment_tables/bam_peak_group_timeseries.tsv"

# bam + group (LeafCutter, MAJIQ, rMATS)
{
    printf 'sample\tbam\tgroup\n'
    for i in 1 2 3; do printf "sample%s\t%s/bams/sample%s.bam\tRef\n" "$i" "$BASE" "$i"; done
    for i in 4 5 6; do printf "sample%s\t%s/bams/sample%s.bam\tAlt\n" "$i" "$BASE" "$i"; done
} > "$BASE/experiment_tables/bam_group.tsv"

# target + control (callpeak_ChIPseq)
{
    printf 'sample\ttarget\tcontrol\n'
    printf "H3K27ac\t%s/chipseq_bams/target1.bam\t%s/chipseq_bams/control1.bam\n" "$BASE" "$BASE"
} > "$BASE/experiment_tables/chipseq_callpeak.tsv"

# R1 + R2 (STARsolo, cellranger_count, kb-nac)
{
    printf 'sample\tR1\tR2\n'
    printf "sample1\t%s/fastqs/sample1_1.fastq.gz\t%s/fastqs/sample1_2.fastq.gz\n" "$BASE" "$BASE"
    printf "sample2\t%s/fastqs/sample2_1.fastq.gz\t%s/fastqs/sample2_2.fastq.gz\n" "$BASE" "$BASE"
} > "$BASE/experiment_tables/R1_R2.tsv"

# fastq + group (Whippet, SUPPA2_diffSplice)
{
    printf 'sample\tfastq\tgroup\n'
    for i in 1 2 3; do printf "sample%s\t%s/fastqs/sample%s_1.fastq.gz,%s/fastqs/sample%s_2.fastq.gz\tRef\n" "$i" "$BASE" "$i" "$BASE" "$i"; done
    for i in 4 5 6; do printf "sample%s\t%s/fastqs/sample%s_1.fastq.gz,%s/fastqs/sample%s_2.fastq.gz\tAlt\n" "$i" "$BASE" "$i" "$BASE" "$i"; done
} > "$BASE/experiment_tables/fastq_group.tsv"

###############################################################################
# bam_list for bam2cram
###############################################################################

printf 'sample1.bam\nsample2.bam\n' > "$BASE/bam_list.txt"

###############################################################################
# Workdirs with source files
###############################################################################

# preprocessing_RNAseq (paired)
mkdir -p "$BASE/workdirs/rnaseq_paired/fastq"
touch "$BASE/workdirs/rnaseq_paired/fastq/sample1_1.fastq.gz"
touch "$BASE/workdirs/rnaseq_paired/fastq/sample1_2.fastq.gz"

# preprocessing_RNAseq (single)
mkdir -p "$BASE/workdirs/rnaseq_single/fastq"
touch "$BASE/workdirs/rnaseq_single/fastq/sample1.fastq.gz"

# preprocessing_RNAseq_long
mkdir -p "$BASE/workdirs/rnaseq_long/fastq"
touch "$BASE/workdirs/rnaseq_long/fastq/sample1.fastq.gz"

# preprocessing_ChIPseq (paired)
mkdir -p "$BASE/workdirs/chipseq_paired/fastq"
touch "$BASE/workdirs/chipseq_paired/fastq/sample1_1.fastq.gz"
touch "$BASE/workdirs/chipseq_paired/fastq/sample1_2.fastq.gz"

# preprocessing_ChIPseq (single)
mkdir -p "$BASE/workdirs/chipseq_single/fastq"
touch "$BASE/workdirs/chipseq_single/fastq/sample1.fastq.gz"

# preprocessing_CLIPseq (HITSCLIP + iCLIPseq)
for d in clipseq_hitsclip clipseq_iclipseq; do
    mkdir -p "$BASE/workdirs/$d/fastq"
    touch "$BASE/workdirs/$d/fastq/sample1.fastq.gz"
done

# callpeak_ATACseq
mkdir -p "$BASE/workdirs/callpeak_atacseq/bowtie2"
touch "$BASE/workdirs/callpeak_atacseq/bowtie2/sample1.sort.rmdup.bam"

# callpeak_ChIPseq
mkdir -p "$BASE/workdirs/callpeak_chipseq"

# differential_ATACseq
mkdir -p "$BASE/workdirs/differential_atacseq"

# footprinting_ATACseq
mkdir -p "$BASE/workdirs/footprinting_atacseq"

# footprinting_timeseries_ATACseq
mkdir -p "$BASE/workdirs/footprinting_ts_atacseq"

# bam2cram (no workdir: directive; uses --directory)
mkdir -p "$BASE/workdirs/bam2cram"
touch "$BASE/workdirs/bam2cram/sample1.bam" "$BASE/workdirs/bam2cram/sample2.bam"

# STARsolo
mkdir -p "$BASE/workdirs/starsolo"

# cellranger_count
mkdir -p "$BASE/workdirs/cellranger/fastq"

# kb-nac
mkdir -p "$BASE/workdirs/kb_nac"

# ONT_plasmid_assembly
mkdir -p "$BASE/workdirs/ont_assembly"
mkdir -p "$BASE/ont_fastq_dirs/sample1"
touch "$BASE/ont_fastq_dirs/sample1/reads.fastq.gz"

# Whippet
mkdir -p "$BASE/workdirs/whippet"

# SUPPA2
mkdir -p "$BASE/workdirs/suppa2"

# LeafCutter
mkdir -p "$BASE/workdirs/leafcutter"

# rMATS
mkdir -p "$BASE/workdirs/rmats"

# MAJIQ
mkdir -p "$BASE/workdirs/majiq"

###############################################################################
# Config YAML files
###############################################################################

mkdir -p "$BASE/configs"

# 1. preprocessing_RNAseq (paired)
cat > "$BASE/configs/config_preprocessing_RNAseq_paired.yaml" << YAML
workdir: ${BASE}/workdirs/rnaseq_paired
samples: ["sample1"]
star_index: ${BASE}/star_index
gtf: ${BASE}/dummy.gtf
layout: "paired"
YAML

# 2. preprocessing_RNAseq (single)
cat > "$BASE/configs/config_preprocessing_RNAseq_single.yaml" << YAML
workdir: ${BASE}/workdirs/rnaseq_single
samples: ["sample1"]
star_index: ${BASE}/star_index
gtf: ${BASE}/dummy.gtf
layout: "single"
YAML

# 3. preprocessing_RNAseq_long
cat > "$BASE/configs/config_preprocessing_RNAseq_long.yaml" << YAML
workdir: ${BASE}/workdirs/rnaseq_long
samples: ["sample1"]
genome_fasta: ${BASE}/dummy.fa
gtf: ${BASE}/dummy.gtf
preset: "splice"
minimap2_options: ""
data_type: "nanopore"
isoquant_options: "--complete_genedb"
YAML

# 4. preprocessing_ChIPseq (paired)
cat > "$BASE/configs/config_preprocessing_ChIPseq_paired.yaml" << YAML
workdir: ${BASE}/workdirs/chipseq_paired
samples: ["sample1"]
bowtie2_index: "${BASE}/bowtie2_index/genome"
layout: "paired"
bowtie2_args: ""
YAML

# 5. preprocessing_ChIPseq (single)
cat > "$BASE/configs/config_preprocessing_ChIPseq_single.yaml" << YAML
workdir: ${BASE}/workdirs/chipseq_single
samples: ["sample1"]
bowtie2_index: "${BASE}/bowtie2_index/genome"
layout: "single"
bowtie2_args: ""
YAML

# 6. preprocessing_CLIPseq (HITSCLIP)
cat > "$BASE/configs/config_preprocessing_CLIPseq_HITSCLIP.yaml" << YAML
workdir: ${BASE}/workdirs/clipseq_hitsclip
samples: ["sample1"]
star_index: ${BASE}/star_index
fastp_args: "-l 20 -3 --trim_front1 5"
outFilterMultimapNmax: 100
normalize_bigwig: false
YAML

# 7. preprocessing_CLIPseq (iCLIPseq)
cat > "$BASE/configs/config_preprocessing_CLIPseq_iCLIPseq.yaml" << YAML
workdir: ${BASE}/workdirs/clipseq_iclipseq
samples: ["sample1"]
star_index: ${BASE}/star_index
fastp_args: "--trim_front1 9 --max_len1 35 --length_required 25"
outFilterMultimapNmax: 1
normalize_bigwig: true
YAML

# 8. callpeak_ATACseq
cat > "$BASE/configs/config_callpeak_ATACseq.yaml" << YAML
workdir: ${BASE}/workdirs/callpeak_atacseq
samples: ["sample1"]
broad: false
broad_cutoff:
genomesize: "mm"
fileformat: "BAM"
qvalue: 0.05
YAML

# 9. callpeak_ChIPseq
cat > "$BASE/configs/config_callpeak_ChIPseq.yaml" << YAML
general:
  workdir: ${BASE}/workdirs/callpeak_chipseq
  experiment_table: ${BASE}/experiment_tables/chipseq_callpeak.tsv
macs2:
  broad: false
  broad_cutoff:
  genomesize: "mm"
  fileformat: "BAM"
  qvalue: 0.05
bedgraphtobigwig:
  assembly: "mm10"
YAML

# 10. differential_ATACseq
cat > "$BASE/configs/config_differential_ATACseq.yaml" << YAML
workdir: ${BASE}/workdirs/differential_atacseq
experiment_table: ${BASE}/experiment_tables/bam_peak_group.tsv
YAML

# 11. footprinting_ATACseq
cat > "$BASE/configs/config_footprinting_ATACseq.yaml" << YAML
workdir: ${BASE}/workdirs/footprinting_atacseq
experiment_table: ${BASE}/experiment_tables/bam_peak_group.tsv
genome_fasta: ${BASE}/dummy.fa
cluster_motifs: ${BASE}/motifs.jaspar
YAML

# 12. footprinting_timeseries_ATACseq
cat > "$BASE/configs/config_footprinting_timeseries_ATACseq.yaml" << YAML
workdir: ${BASE}/workdirs/footprinting_ts_atacseq
experiment_table: ${BASE}/experiment_tables/bam_peak_group_timeseries.tsv
genome_fasta: ${BASE}/dummy.fa
cluster_motifs: ${BASE}/motifs.jaspar
YAML

# 13. bam2cram
cat > "$BASE/configs/config_bam2cram.yaml" << YAML
reference: ${BASE}/dummy.fa
bam_list: ${BASE}/bam_list.txt
YAML

# 14. STARsolo
cat > "$BASE/configs/config_STARsolo.yaml" << YAML
workdir: ${BASE}/workdirs/starsolo
experiment_table: ${BASE}/experiment_tables/R1_R2.tsv
star_index: ${BASE}/star_index
barcode_whitelist: ${BASE}/barcode_whitelist.txt
soloUMIlen: 12
soloCellFilter: EmptyDrops_CR
clipAdapterType: CellRanger4
outFilterScoreMin: 30
soloCBmatchWLtype: 1MM_multi_Nbase_pseudocounts
soloUMIfiltering: MultiGeneUMI_CR
soloUMIdedup: 1MM_CR
soloBarcodeReadLength: 150
YAML

# 15. cellranger_count
cat > "$BASE/configs/config_cellranger_count.yaml" << YAML
workdir: ${BASE}/workdirs/cellranger
experiment_table: ${BASE}/experiment_tables/R1_R2.tsv
transcriptome: ${BASE}/transcriptome
create_bam: "false"
YAML

# 16. kb-nac
cat > "$BASE/configs/config_kb_nac.yaml" << YAML
workdir: ${BASE}/workdirs/kb_nac
experiment_table: ${BASE}/experiment_tables/R1_R2.tsv
dna_fasta: ${BASE}/dummy.fa
gtf: ${BASE}/dummy.gtf
technology: 10xv3
YAML

# 17. ONT_plasmid_assembly
cat > "$BASE/configs/config_ONT_plasmid_assembly.yaml" << YAML
workdir: ${BASE}/workdirs/ont_assembly
samples:
  sample1: "${BASE}/ont_fastq_dirs/sample1"
reference: ${BASE}/dummy.fa
YAML

# 18. Whippet
cat > "$BASE/configs/config_Whippet.yaml" << YAML
workdir: ${BASE}/workdirs/whippet
experiment_table: ${BASE}/experiment_tables/fastq_group.tsv
whippet_index: ${BASE}/whippet_index.jls
YAML

# 19. SUPPA2_diffSplice
cat > "$BASE/configs/config_SUPPA2_diffSplice.yaml" << YAML
workdir: ${BASE}/workdirs/suppa2
experiment_table: ${BASE}/experiment_tables/fastq_group.tsv
salmon_index: ${BASE}/salmon_index
gtf: ${BASE}/dummy.gtf
YAML

# 20. LeafCutter
cat > "$BASE/configs/config_LeafCutter.yaml" << YAML
workdir: ${BASE}/workdirs/leafcutter
experiment_table: ${BASE}/experiment_tables/bam_group.tsv
gtf: ${BASE}/dummy.gtf
minimum_anchor_length: 6
minimum_intron_length: 70
maximum_intron_length: 500000
strand: XS
minimum_reads: 10
min_coverage: 10
min_samples_per_intron: 1
min_samples_per_group: 3
FDR: 0.05
YAML

# 21. rMATS
cat > "$BASE/configs/config_rMATS.yaml" << YAML
workdir: ${BASE}/workdirs/rmats
experiment_table: ${BASE}/experiment_tables/bam_group.tsv
gtf: ${BASE}/dummy.gtf
layout: paired
readLength: 101
cstat: 0.1
anchorLength: 6
mil: 0
mel: 10000
YAML

# 22. MAJIQ
cat > "$BASE/configs/config_MAJIQ.yaml" << YAML
container: dummy_container.sif
workdir: ${BASE}/workdirs/majiq
experiment_table: ${BASE}/experiment_tables/bam_group.tsv
gff: ${BASE}/dummy.gff3
genome: hg19
strandness: None
YAML

echo "Test data generated successfully in $BASE"
