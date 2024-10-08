
'''
Snakefile for LeafCutter

Usage:
    snakemake -s LeafCutter.smk --configfile config.yaml --cores <int> --use-singularity
'''

'''experiment_table
sample	bam	group
SRR536342	/rhome/naotok/bigdata/Shiba/GEO/SRP014759/star/SRR536342/SRR536342_Aligned.out.bam	Ref
SRR536344	/rhome/naotok/bigdata/Shiba/GEO/SRP014759/star/SRR536344/SRR536344_Aligned.out.bam	Ref
SRR536346	/rhome/naotok/bigdata/Shiba/GEO/SRP014759/star/SRR536346/SRR536346_Aligned.out.bam	Ref
SRR536348	/rhome/naotok/bigdata/Shiba/GEO/SRP014759/star/SRR536348/SRR536348_Aligned.out.bam	Alt
SRR536350	/rhome/naotok/bigdata/Shiba/GEO/SRP014759/star/SRR536350/SRR536350_Aligned.out.bam	Alt
SRR536352	/rhome/naotok/bigdata/Shiba/GEO/SRP014759/star/SRR536352/SRR536352_Aligned.out.bam	Alt
'''

workdir: config['workdir']

import os
import pandas as pd

def get_bam_files(experiment_table):
    df = pd.read_csv(experiment_table, sep = '\t')
    return {sample: bam for sample, bam in zip(df['sample'], df['bam'])}

def get_groups(experiment_table):
    df = pd.read_csv(experiment_table, sep = '\t')
    return {sample: group for sample, group in zip(df['sample'], df['group'])}

sample_bam_dict = get_bam_files(config['experiment_table'])
sample_group_dict = get_groups(config['experiment_table'])
samples = sample_bam_dict.keys()

rule all:
    input:
        cluster_significance = "differential_cluster_significance.txt",
        effect_sizes = "differential_effect_sizes.txt",
        leafviz_rdata = "leafviz.Rdata",
        cluster_classification = "cluster_classifications.tsv",
        cluster_summary = "cluster_summary.tsv",

rule makeJuncfile:
    output:
        juncfile = "juncfiles.txt"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/makeJuncfiles.txt"
    log:
        "log/leafcutter/makeJuncfiles.log"
    run:
        with open(output.juncfile, 'w') as f:
            for sample in samples:
                f.write(f"junc/{sample}.junc\n")

rule makeGroupfile:
    output:
        groupfile = "groupfile.txt"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/makeGroupfile.txt"
    log:
        "log/leafcutter/makeGroupfile.log"
    run:
        with open(output.groupfile, 'w') as f:
            for sample in samples:
                f.write(f"{sample}\t{sample_group_dict[sample]}\n")

rule regtools:
    wildcard_constraints:
        sample = "|".join(samples)
    container:
        "docker://griffithlab/regtools:release-1.0.0"
    input:
        bam = lambda wildcards: sample_bam_dict[wildcards.sample],
    output:
        junc = "junc/{sample}.junc"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/Regtools_{sample}.txt"
    log:
        "log/leafcutter/Regtools_{sample}.log"
    shell:
        """
        mkdir -p junc/ && \
        regtools junctions extract \
        -s {config[strand]} \
        -a {config[minimum_anchor_length]} \
        -m {config[minimum_intron_length]} \
        -M {config[maximum_intron_length]} \
        {input.bam} \
        -o {output.junc} \
        >& {log}
        """

rule intronClustering:
    container:
        "docker://naotokubota/leafcutter:0.2.9"
    input:
        junc = expand("junc/{sample}.junc", sample = samples),
        juncfile = "juncfiles.txt"
    output:
        clusters = "clusters_perind_numers.counts.gz"
    params:
        output_prefix = "clusters"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/intronClustering.txt"
    log:
        "log/leafcutter/intronClustering.log"
    shell:
        """
        python /opt/leafcutter/leafcutter/scripts/leafcutter_cluster_regtools.py \
        -j {input.juncfile} \
        -m {config[minimum_reads]} \
        -l {config[maximum_intron_length]} \
        -o {params.output_prefix} \
        >& {log}
        """

rule gtf2exon:
    container:
        "docker://naotokubota/leafcutter:0.2.9"
    output:
        exon = "exons.txt.gz"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/gtf2exon.txt"
    log:
        "log/leafcutter/gtf2exon.log"
    shell:
        """
        /opt/leafcutter/leafcutter/scripts/gtf_to_exons.R \
        {config[gtf]} \
        {output.exon} \
        >& {log}
        """

rule differentialAnalysis:
    container:
        "docker://naotokubota/leafcutter:0.2.9"
    input:
        clusters = "clusters_perind_numers.counts.gz",
        groupfile = "groupfile.txt",
        exon = "exons.txt.gz"
    output:
        cluster_significance = "differential_cluster_significance.txt",
        effect_sizes = "differential_effect_sizes.txt"
    params:
        output_prefix = "differential"
    threads:
        workflow.cores
    benchmark:
        "benchmark/leafcutter/differentialAnalysis.txt"
    log:
        "log/leafcutter/differentialAnalysis.log"
    shell:
        """
        /opt/leafcutter/leafcutter/scripts/leafcutter_ds.R \
        -p {threads} \
        -c {config[min_coverage]} \
        -i {config[min_samples_per_intron]} \
        -e {input.exon} \
        -o {params.output_prefix} \
        -g {config[min_samples_per_group]} \
        {input.clusters} \
        {input.groupfile} \
        >& {log}
        """

rule makeAnnotationCodes:
    container:
        "docker://naotokubota/leafcutter:0.2.9"
    output:
        all_exons = "annotation_codes/annotation_codes_all_exons.txt.gz",
        all_introns = "annotation_codes/annotation_codes_all_introns.bed.gz",
        fiveprime = "annotation_codes/annotation_codes_fiveprime.bed.gz",
        threeprime = "annotation_codes/annotation_codes_threeprime.bed.gz"
    params:
        output_prefix = "annotation_codes/annotation_codes"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/makeAnnotationCodes.txt"
    log:
        "log/leafcutter/makeAnnotationCodes.log"
    shell:
        """
        mkdir -p annotation_codes/ && \
        /opt/leafcutter/leafcutter/leafviz/gtf2leafcutter.pl \
        -o {params.output_prefix} \
        {config[gtf]} \
        >& {log}
        """

rule prepareResults:
    container:
        "docker://naotokubota/leafcutter:0.2.9"
    input:
        clusters = "clusters_perind_numers.counts.gz",
        groupfile = "groupfile.txt",
        cluster_significance = "differential_cluster_significance.txt",
        effect_sizes = "differential_effect_sizes.txt",
        all_exons = "annotation_codes/annotation_codes_all_exons.txt.gz",
        all_introns = "annotation_codes/annotation_codes_all_introns.bed.gz",
        fiveprime = "annotation_codes/annotation_codes_fiveprime.bed.gz",
        threeprime = "annotation_codes/annotation_codes_threeprime.bed.gz"
    params:
        annotation_codes_prefix = "annotation_codes/annotation_codes"
    output:
        leafviz_rdata = "leafviz.Rdata"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/prepareResults.txt"
    log:
        "log/leafcutter/prepareResults.log"
    shell:
        """
        /opt/leafcutter/leafcutter/leafviz/prepare_results.R \
        -m {input.groupfile} \
        {input.clusters} \
        {input.cluster_significance} \
        {input.effect_sizes} \
        {params.annotation_codes_prefix} \
        -o {output.leafviz_rdata} \
        >& {log}
        """

rule classifyClusters:
    container:
        "docker://naotokubota/leafcutter:0.2.9"
    input:
        leafviz_rdata = "leafviz.Rdata"
    output:
        cluster_classification = "cluster_classifications.tsv",
        cluster_summary = "cluster_summary.tsv"
    params:
        output_prefix = "cluster"
    threads:
        1
    benchmark:
        "benchmark/leafcutter/classifyClusters.txt"
    log:
        "log/leafcutter/classifyClusters.log"
    shell:
        """
        Rscript /opt/leafcutter/leafcutter/leafcutter/R/classify_clusters.R \
        -o {params.output_prefix} \
        {input.leafviz_rdata} \
        >& {log}
        """
