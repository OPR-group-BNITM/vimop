#!/usr/bin/env bash

nextflow main.nf \
    --fastq /Users/nils.petersen/dev/metagenomics_pipeline/nanoflow_analysis_data/TEST-DATA/miniset/pass \
    --sample_sheet samplesheet.csv 
