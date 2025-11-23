# Quality control: FastQC on raw FASTQ files

rule fastqc_raw:
    """
    Run FastQC on raw FASTQ files (R1 and R2) for each sample.
    """
    input:
        r1=lambda wc: fq(wc.sample, 1),
        r2=lambda wc: fq(wc.sample, 2)
    output:
        html_r1=f"{QC_DIR}/fastqc/{{sample}}_R1_fastqc.html",
        zip_r1=f"{QC_DIR}/fastqc/{{sample}}_R1_fastqc.zip",
        html_r2=f"{QC_DIR}/fastqc/{{sample}}_R2_fastqc.html",
        zip_r2=f"{QC_DIR}/fastqc/{{sample}}_R2_fastqc.zip",
    threads:
        config["threads"]["fastqc"]
    conda:
        "envs/germline.yaml"
    shell:
        r"""
        mkdir -p {QC_DIR}/fastqc

        fastqc -t {threads} \
            -o {QC_DIR}/fastqc \
            {input.r1} {input.r2}
        """
