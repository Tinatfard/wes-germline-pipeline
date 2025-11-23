# Mapping: BWA-MEM + sorted BAM + index

rule bwa_mem_sort_index:
    """
    Map reads with BWA-MEM, sort BAM, and index.
    """
    input:
        ref = REF,
        r1 = lambda wc: fq(wc.sample, 1),
        r2 = lambda wc: fq(wc.sample, 2)
    output:
        bam = f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai = f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai"
    threads:
        config["threads"]["bwa"]
    conda:
        "envs/germline.yaml"
    shell:
        r"""
        mkdir -p {OUTDIR}/bam

        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} \
            | samtools sort -o {output.bam} -

        samtools index {output.bam}
        """
