# Variant calling: GATK HaplotypeCaller (per-sample GVCF) + GenotypeGVCFs (cohort VCF)

# Uses global variables from Snakefile:
# REF, OUTDIR, VAR_DIR, SAMPLES, config


rule haplotypecaller:
    """
    Call germline variants per-sample using GATK HaplotypeCaller in GVCF mode.
    Produces one .g.vcf.gz + .tbi per sample.
    """
    input:
        bam = f"{OUTDIR}/bam/{{sample}}.sorted.bam",
        bai = f"{OUTDIR}/bam/{{sample}}.sorted.bam.bai",
        ref = REF
    output:
        gvcf = f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz",
        gvcf_idx = f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz.tbi"
    threads:
        config["threads"]["gatk"]
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


rule genotypegvcfs:
    """
    Jointly genotype all per-sample GVCFs into a single cohort VCF.
    """
    input:
        gvcfs = expand(f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz", sample=SAMPLES),
        gvcfs_idx = expand(f"{OUTDIR}/gvcf/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES),
        ref = REF
    output:
        vcf = f"{VAR_DIR}/cohort.vcf.gz",
        vcf_idx = f"{VAR_DIR}/cohort.vcf.gz.tbi"
    threads:
        config["threads"]["gatk"]
    params:
        extra = config.get("gatk", {}).get("genotypegvcfs_extra", "")
    conda:
        "envs/germline.yaml"
    shell:
        r"""
        mkdir -p {VAR_DIR}

        gatk GenotypeGVCFs \
            -R {input.ref} \
            { " ".join(f"-V {g}" for g in input.gvcfs) } \
            -O {output.vcf} \
            {params.extra}

        tabix -p vcf {output.vcf}
        """

