'''
Snakefile for de novo transcriptome assembly from short-read RNA-seq data (paired-end and single-end)

Usage:
    snakemake -s denovo_assembly_RNAseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity

config.yaml must contain:
    layout: "paired" or "single"
    max_memory: e.g. "32G"
    swissprot_db: path to Swiss-Prot BLAST database
'''

## DON'T CHANGE BELOW THIS LINE ##

import re

import os

workdir: config["workdir"]
samples = config["samples"]
layout = config["layout"]
swissprot_db = config["swissprot_db"]
swissprot_db_dir = os.path.dirname(os.path.abspath(swissprot_db))
max_memory = config["max_memory"]

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

include: "common/containers.smk"

if layout == "paired":

    rule all:
        input:
            trinity = "trinity/Trinity.fasta",
            trinity_stats = "trinity/Trinity_stats.txt",
            transdecoder_pep = "transdecoder/Trinity.fasta.transdecoder.pep",
            transdecoder_bed = "transdecoder/Trinity.fasta.transdecoder.bed",
            transdecoder_gff3 = "transdecoder/Trinity.fasta.transdecoder.gff3",
            blastx = "blast/blastx_swissprot.outfmt6",
            annotation = "annotation/annotation_summary.tsv"

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            R1 = "fastq/{sample}_1.fastq.gz",
            R2 = "fastq/{sample}_2.fastq.gz"
        output:
            R1 = "fastp/{sample}_1.fastq.gz",
            R2 = "fastp/{sample}_2.fastq.gz",
            json = "fastp/log/{sample}.json"
        params:
            html = "fastp/log/{sample}.html"
        threads:
            8
        benchmark:
            "benchmark/fastp_{sample}.txt"
        log:
            "log/fastp_{sample}.log"
        shell:
            "fastp -i {input.R1} -I {input.R2} "
            "-o {output.R1} -O {output.R2} -w {threads} "
            "-h {params.html} -j {output.json} >& {log} && "
            "for f in {output.R1} {output.R2}; do "
            "  zcat $f | awk 'NR%4==1{{print $1; next}}{{print}}' | gzip > $f.tmp && mv $f.tmp $f; "
            "done"

    rule trinity:
        container:
            CONTAINERS["trinity"]
        input:
            R1 = expand("fastp/{sample}_1.fastq.gz", sample = samples),
            R2 = expand("fastp/{sample}_2.fastq.gz", sample = samples)
        output:
            fasta = "trinity/Trinity.fasta"
        params:
            left = lambda wildcards, input: ",".join(input.R1),
            right = lambda wildcards, input: ",".join(input.R2),
            outdir = "trinity"
        threads:
            workflow.cores
        benchmark:
            "benchmark/trinity.txt"
        log:
            "log/trinity.log"
        shell:
            "Trinity --seqType fq "
            "--left {params.left} --right {params.right} "
            "--max_memory {max_memory} --CPU {threads} "
            "--output {params.outdir} >& {log} && "
            "mv {params.outdir}.Trinity.fasta {output.fasta} && "
            "mv {params.outdir}.Trinity.fasta.gene_trans_map {params.outdir}/Trinity.fasta.gene_trans_map"

else:

    rule all:
        input:
            trinity = "trinity/Trinity.fasta",
            trinity_stats = "trinity/Trinity_stats.txt",
            transdecoder_pep = "transdecoder/Trinity.fasta.transdecoder.pep",
            transdecoder_bed = "transdecoder/Trinity.fasta.transdecoder.bed",
            transdecoder_gff3 = "transdecoder/Trinity.fasta.transdecoder.gff3",
            blastx = "blast/blastx_swissprot.outfmt6",
            annotation = "annotation/annotation_summary.tsv"

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            R1 = "fastq/{sample}.fastq.gz"
        output:
            R1 = "fastp/{sample}.fastq.gz",
            json = "fastp/log/{sample}.json"
        params:
            html = "fastp/log/{sample}.html"
        threads:
            8
        benchmark:
            "benchmark/fastp_{sample}.txt"
        log:
            "log/fastp_{sample}.log"
        shell:
            "fastp -i {input.R1} "
            "-o {output.R1} -w {threads} "
            "-h {params.html} -j {output.json} >& {log}"

    rule trinity:
        container:
            CONTAINERS["trinity"]
        input:
            R1 = expand("fastp/{sample}.fastq.gz", sample = samples)
        output:
            fasta = "trinity/Trinity.fasta"
        params:
            single = lambda wildcards, input: ",".join(input.R1),
            outdir = "trinity"
        threads:
            workflow.cores
        benchmark:
            "benchmark/trinity.txt"
        log:
            "log/trinity.log"
        shell:
            "Trinity --seqType fq "
            "--single {params.single} "
            "--max_memory {max_memory} --CPU {threads} "
            "--output {params.outdir} >& {log} && "
            "mv {params.outdir}.Trinity.fasta {output.fasta} && "
            "mv {params.outdir}.Trinity.fasta.gene_trans_map {params.outdir}/Trinity.fasta.gene_trans_map"

rule assembly_stats:
    container:
        CONTAINERS["trinity"]
    input:
        fasta = "trinity/Trinity.fasta"
    output:
        stats = "trinity/Trinity_stats.txt"
    threads:
        1
    benchmark:
        "benchmark/trinity_stats.txt"
    log:
        "log/trinity_stats.log"
    shell:
        "/usr/local/bin/util/TrinityStats.pl {input.fasta} > {output.stats} 2> {log}"

rule transdecoder_longorfs:
    container:
        CONTAINERS["transdecoder"]
    input:
        fasta = "trinity/Trinity.fasta"
    output:
        longorfs = "transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep"
    params:
        outdir = "transdecoder"
    threads:
        1
    benchmark:
        "benchmark/transdecoder_longorfs.txt"
    log:
        "log/transdecoder_longorfs.log"
    shell:
        "cd {params.outdir} && "
        "TransDecoder.LongOrfs -t ../trinity/Trinity.fasta >& ../{log}"

rule blastx_swissprot:
    container:
        CONTAINERS["blast"]
    input:
        fasta = "trinity/Trinity.fasta"
    output:
        outfmt6 = "blast/blastx_swissprot.outfmt6"
    params:
        db = swissprot_db
    threads:
        workflow.cores
    benchmark:
        "benchmark/blastx_swissprot.txt"
    log:
        "log/blastx_swissprot.log"
    shell:
        "blastx -query {input.fasta} -db {params.db} "
        "-max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' "
        "-evalue 1e-5 -num_threads {threads} "
        "-out {output.outfmt6} >& {log}"

rule blastp_swissprot:
    container:
        CONTAINERS["blast"]
    input:
        pep = "transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        outfmt6 = "blast/blastp_swissprot.outfmt6"
    params:
        db = swissprot_db
    threads:
        workflow.cores
    benchmark:
        "benchmark/blastp_swissprot.txt"
    log:
        "log/blastp_swissprot.log"
    shell:
        "blastp -query {input.pep} -db {params.db} "
        "-max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' "
        "-evalue 1e-5 -num_threads {threads} "
        "-out {output.outfmt6} >& {log}"

rule transdecoder_predict:
    container:
        CONTAINERS["transdecoder"]
    input:
        fasta = "trinity/Trinity.fasta",
        longorfs = "transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep",
        blastp = "blast/blastp_swissprot.outfmt6"
    output:
        pep = "transdecoder/Trinity.fasta.transdecoder.pep",
        bed = "transdecoder/Trinity.fasta.transdecoder.bed",
        gff3 = "transdecoder/Trinity.fasta.transdecoder.gff3"
    params:
        outdir = "transdecoder"
    threads:
        1
    benchmark:
        "benchmark/transdecoder_predict.txt"
    log:
        "log/transdecoder_predict.log"
    shell:
        "cd {params.outdir} && "
        "TransDecoder.Predict -t ../trinity/Trinity.fasta "
        "--retain_blastp_hits ../blast/blastp_swissprot.outfmt6 "
        "--single_best_only >& ../{log}"

rule annotation_summary:
    input:
        blastx = "blast/blastx_swissprot.outfmt6",
        transdecoder_pep = "transdecoder/Trinity.fasta.transdecoder.pep"
    output:
        summary = "annotation/annotation_summary.tsv"
    threads:
        1
    benchmark:
        "benchmark/annotation_summary.txt"
    log:
        "log/annotation_summary.log"
    shell:
        """
        (
        echo -e "transcript_id\\tgene_id\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tdescription" > {output.summary} && \
        cat {input.blastx} >> {output.summary} && \
        echo "" >> {output.summary} && \
        echo "# TransDecoder predicted ORFs (from .pep headers):" >> {output.summary} && \
        grep "^>" {input.transdecoder_pep} | sed 's/^>//' >> {output.summary}
        ) >& {log}
        """

## DON'T CHANGE ABOVE THIS LINE ##
