rule bwa_mem:
    """
    Align reads with BWA-MEM and write an unsorted BAM.
    """
    input:
        r1=lambda wc: fq(wc.sample, 1),
        r2=lambda wc: fq(wc.sample, 2),
        ref=REF
    output:
        bam=temp(f"{OUTDIR}/bam/{{sample}}.unsorted.bam")
    threads: config["threads"]["bwa"]
    conda:
        "../envs/germline.yaml"
    shell:
        r"""
        mkdir -p {OUTDIR}/bam
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
          samtools view -bS - > {output.bam}
        """

rule sort_bam:
    """
    Sort BAM and index it.
    """
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.unsorted.bam"
    output:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai"
    threads: 2
    conda:
        "../envs/germline.yaml"
    shell:
        r"""
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index {output.bam}
        """
