rule haplotypecaller:
    """
    Call germline variants using GATK HaplotypeCaller in GVCF mode.
    Produces one .g.vcf.gz per sample.
    """
    input:
        bam = f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai = f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai",
        ref = REF
    output:
        gvcf = f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz",
        gvcf_idx = f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz.tbi"
    threads:
        config["threads"].get("gatk", 4)
    params:
        extra = config.get("gatk", {}).get("haplotypecaller_extra", "")
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

        tabix -p vcf {output.gvcf}
        """

