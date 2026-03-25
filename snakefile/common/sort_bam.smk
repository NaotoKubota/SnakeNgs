rule sort:
    container:
        CONTAINERS["samtools"]
    input:
        "star/{sample}/{sample}_Aligned.out.sam"
    output:
        "star/{sample}/{sample}_Aligned.out.bam"
    params:
        index_option = "-c" if check_chromosome_length(star_index) else "-b"
    threads:
        workflow.cores / 4
    benchmark:
        "benchmark/samtools_{sample}.txt"
    log:
        "log/samtools_{sample}.log"
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input} >& {log} && "
        "samtools index {params.index_option} {output} && "
        "rm -rf {input}"
