#!/usr/bin/env bash

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.


# Directory of this script
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
datadir=$script_dir/data/custom_db_test

cd "$script_dir/.."

nextflow main.nf \
    --custom_db_do_build \
    --custom_db_outpath custom_db_testout \
    --custom_db_contaminants_version 1.0.0 \
    --custom_db_contaminants_description "This is a dummy contaminants data base to test the custom data base module." \
    --custom_db_contaminants_input_path "${datadir}/contaminants" \
    --custom_db_centrifuge_version 1.0.0 \
    --custom_db_centrifuge_description "This is a dummy centrifuge data base to test the custom data base module." \
    --custom_db_centrifuge_fasta "${datadir}/centrifuge/input.fasta" \
    --custom_db_virus_yaml "${datadir}/virus/testset.yaml" \
    --custom_db_virus_fasta "${datadir}/virus/ALL.fasta" \
    --custom_db_centrifuge_version 1.0.0 \
    --custom_db_centrifuge_description "This is a dummy virus reference data base to test the custom data base module." \
    --custom_db_min_disk_space_work_gb 50 \
    --custom_db_min_disk_space_out_gb 50
