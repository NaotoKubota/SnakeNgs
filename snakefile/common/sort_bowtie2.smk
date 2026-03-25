rule sort:
    container:
        CONTAINERS["samtools"]
    input:
        "bowtie2/{sample}.sam"
    output:
        temp("bowtie2/{sample}.sort.bam")
    threads:
        8
    benchmark:
        "benchmark/samtools_{sample}.txt"
    log:
        "log/samtools_{sample}.log"
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input} >& {log} && "
        "samtools index {output}"
