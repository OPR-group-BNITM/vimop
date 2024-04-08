#!/usr/bin/env bash

 nextflow main.nf \
    --fastq /Users/nils.petersen/dev/metagenomics_pipeline/nanoflow_analysis_data/TEST-DATA/miniset/pass \
    --sample_sheet samplesheet.csv 

# nextflow run main.nf \
#     --fastq /Volumes/Samsung_T5/assembly_testset/HAM-2021-RUN013/pass \
#     --sample_sheet /Volumes/Samsung_T5/assembly_testset/HAM-2021-RUN013/pass/samplesheet.csv \
#     -work-dir /Volumes/Samsung_T5/metagenomics_work \
#     --out_dir /Volumes/Samsung_T5/metagenomics_out

