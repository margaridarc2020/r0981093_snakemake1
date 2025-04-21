#!/usr/bin/env bash
export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2025/bin/:$PATH
# symlink fastq files to 000.fastq

mkdir 000.fastq
ln -sf /staging/leuven/stg_00079/teaching/data_manual_snpcall/*.fastq 000.fastq/

# Run Snakemake with your snakefile
snakemake -c 24 -s snakefile