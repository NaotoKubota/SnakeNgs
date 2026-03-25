rule makeRefFlat:
    container:
        CONTAINERS["ucsc_gtftogenepred"]
    output:
        refFlat = "metrics/refFlat.txt"
    threads:
        1
    benchmark:
        "benchmark/refFlat.txt"
    log:
        "log/refFlat.log"
    shell:
        "gtfToGenePred -genePredExt -ignoreGroupsWithoutExons -geneNameAsName2 {gtf} refFlat.tmp >& {log} && "
        "cat refFlat.tmp | awk -F'\t' -v OFS='\t' '$4 != 0{{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}$4 == 0{{print $12,$1,$2,$3,1,$5,$6,$7,$8,$9,$10}}' > {output.refFlat} && "
        "rm -rf refFlat.tmp"

rule makeRibosomalInterval:
    output:
        ribosomalInterval = "metrics/ribosomal_interval.txt"
    threads:
        1
    benchmark:
        "benchmark/ribosomal_interval.txt"
    shell:
        '''
        cat {star_index}/chrNameLength.txt | awk -F"\t" -v OFS="\t" '{{print "@SQ","SN:"$1,"LN:"$2}}' > {output.ribosomalInterval} && \
        cat {gtf} | grep -e 'gene_type "rRNA"' -e 'gene_biotype "rRNA"' | \
        awk -F"\t" -v OFS="\t" '$3 == "transcript" && $4 != 1{{print $1,$4-1,$5,$7,$9}}$3 == "transcript" && $4 == 1{{print $1,$4,$5,$7,$9}}' \
        >> {output.ribosomalInterval}
        '''

rule CollectRnaSeqMetrics:
    container:
        CONTAINERS["picard"]
    input:
        bam = "star/{sample}/{sample}_Aligned.out.bam",
        refFlat = "metrics/refFlat.txt",
        ribosomalInterval = "metrics/ribosomal_interval.txt"
    output:
        RnaSeqMetrics = "metrics/{sample}.picard.analysis.CollectRnaSeqMetrics"
    threads:
        1
    benchmark:
        "benchmark/picard_CollectRnaSeqMetrics_{sample}.txt"
    log:
        "log/picard_CollectRnaSeqMetrics_{sample}.log"
    shell:
        "picard CollectRnaSeqMetrics -I {input.bam} -O {output.RnaSeqMetrics} --REF_FLAT {input.refFlat} --STRAND_SPECIFICITY NONE --RIBOSOMAL_INTERVALS {input.ribosomalInterval} >& {log}"
