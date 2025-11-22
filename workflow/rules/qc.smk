rule fastqc:
    input:
        r1=lambda wc: fq(wc.sample, 1),
        r2=lambda wc: fq(wc.sample, 2)
    output:
        html1=lambda wc: f"{QC_DIR}/fastqc/{wc.sample}_R1_fastqc.html",
        html2=lambda wc: f"{QC_DIR}/fastqc/{wc.sample}_R2_fastqc.html"
    threads: config["threads"]["fastqc"]
    conda:
        "../envs/germline.yaml"
    shell:
        r"""
        mkdir -p {QC_DIR}/fastqc
        fastqc {input.r1} {input.r2} -o {QC_DIR}/fastqc
        """
