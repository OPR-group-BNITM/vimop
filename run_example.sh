#!/usr/bin/env bash

nextflow main.nf --fastq /Volumes/Samsung_T5/nanoflow_testsets/CKY-2023-RUN091/pass

# nextflow main.nf \
#     --fastq /Users/nils.petersen/dev/metagenomics_pipeline/nextflow_testdata/tinytest/pass \
#     --sample_sheet  /Users/nils.petersen/dev/metagenomics_pipeline/nextflow_testdata/tinytest/samplesheet.csv 

#  nextflow main.nf \
#     --fastq /Users/nils.petersen/dev/metagenomics_pipeline/nanoflow_analysis_data/TEST-DATA/miniset/pass \
#     --sample_sheet samplesheet.csv 

# nextflow main.nf \
#     --fastq /Users/nils.petersen/dev/metagenomics_pipeline/nextflow_testdata/HAM-2021-RUN013/pass \
#     --sample_sheet  /Users/nils.petersen/dev/metagenomics_pipeline/nextflow_testdata/HAM-2021-RUN013/samplesheet.csv 

# nextflow main.nf \
#     --fastq /Users/nils.petersen/dev/metagenomics_pipeline/nanoflow_analysis_data/TEST-DATA/CKY-2022-RUN051/pass

# nextflow run main.nf \
#     --fastq /Volumes/Samsung_T5/assembly_testset/HAM-2021-RUN013/pass \
#     --sample_sheet /Volumes/Samsung_T5/assembly_testset/HAM-2021-RUN013/pass/samplesheet.csv \
#     -work-dir /Volumes/Samsung_T5/metagenomics/big_work \
#     --out_dir /Volumes/Samsung_T5/metagenomics/big_out

# nextflow run main.nf \
#     --fastq /Volumes/Samsung_T5/assembly_testset/tinytest/pass \
#     --sample_sheet /Volumes/Samsung_T5/assembly_testset/tinytest/pass/samplesheet.csv \
#     -work-dir /Volumes/Samsung_T5/metagenomics_tiny_work \
#     --out_dir /Volumes/Samsung_T5/metagenomics_tiny_out
