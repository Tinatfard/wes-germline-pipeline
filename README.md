# wes-germline-pipeline
This repository contains a small Snakemake workflow for germline variant calling from whole exome sequencing (WES) data.
The goal is to provide a clear and modular example of a basic variant calling pipeline using commonly used tools.

The workflow includes:

Quality control (FastQC)

Alignment with BWA-MEM

Sorting and indexing

Variant calling with GATK HaplotypeCaller (gVCF mode)

The pipeline is designed to be easy to understand and extend.
