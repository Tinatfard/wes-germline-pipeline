rule fastqc:
    """
    Run FastQC on raw FASTQ files.
    """
    input:
        r1=lambda wc: get_fastq(wc.sample, 1),
        r2=lambda wc: get_fastq(wc.sample, 2)
    output:
        html1=f"{QC_DIR}/fastqc/{{sample}}_R1_fastqc.html",
        html2=f"{QC_DIR}/fastqc/{{sample}}_R2_fastqc.html"
    threads: config["threads"].get("fastqc", 2)
    conda:
        "envs/germline.yaml"
    shell:
        r"""
        mkdir -p {QC_DIR}/fastqc
        fastqc {input.r1} {input.r2} -o {QC_DIR}/fastqc
        """
