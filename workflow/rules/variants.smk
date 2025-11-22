rule haplotypecaller:
    """
    Call germline variants with GATK HaplotypeCaller in GVCF mode.
    """
    input:
        bam=f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai",
        ref=REF
    output:
        gvcf=f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz"
    threads: config["threads"].get("gatk", 4)
    params:
        extra=config.get("gatk", {}).get("haplotypecaller_extra", "")
    conda:
        "envs/germline.yaml"
    shell:
        r"""
        mkdir -p {OUTDIR}/gvcf
        gatk HaplotypeCaller \
          -R {input.ref} \
          -I {input.bam} \
          -O {output.gvcf} \
          -ERC GVCF \
          {params.extra}
        """
