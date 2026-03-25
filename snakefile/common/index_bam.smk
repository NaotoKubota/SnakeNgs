rule index:
    container:
        CONTAINERS["samtools"]
    input:
        "bowtie2/{sample}.sort.rmdup.bam"
    output:
        "bowtie2/{sample}.sort.rmdup.bam.bai"
    shell:
        "samtools index {input}"
