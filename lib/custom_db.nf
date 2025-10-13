// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


nextflow.enable.dsl = 2


process build_contaminants_db {
    label "general"
    cpus 1
    input:
        path("contaminant_dir")
    output:
        path("contaminants")
    """
    set -euo pipefail
    shopt -s nullglob extglob

    # check fasta file validity
    for fname in contaminant_dir/*.{fna,fna.gz,fasta,fasta.gz}
    do
        if ! seqtk seq "\$fname" > /dev/null 2>&1
        then
            echo "Error: invalid FASTA file (\$fname)" >&2
            exit 1
        fi
    done

    # copy files (only if any exist)
    mkdir -p contaminants
    files=( contaminant_dir/*.{fna,fna.gz,fasta,fasta.gz} )
    if [ "\${#files[@]}" -gt 0 ]
    then
        cp "\${files[@]}" contaminants/
    fi

    cd contaminants

    # gzip any uncompressed files
    for fname in *.{fna,fasta}
    do
        gzip -f "\$fname"
    done

    # write YAML without quotes around keys; normalize keys to safe chars
    {
        echo "version: ${params.custom_db_contaminants_version}"
        echo "description: ${params.custom_db_contaminants_description}"
        echo "filters:"
        for fn in *.{fna.gz,fasta.gz}
        do
            filter_name=\${fn%.@(fasta|fna).gz}
            # normalize: keep A–Z, a–z, 0–9, dot, underscore, dash; others -> '_'
            filter_name=\${filter_name//[^A-Za-z0-9._-]/_}
            echo "  \${filter_name}: \${fn}"
        done
    } > contaminants.yaml
    """
}


process build_virus_db {
    label "general"
    cpus 1
    input:
        tuple path("virus_input_config.yaml"), path("ALL.fasta")
    output:
        path("virus")
    """
    set -euo pipefail

    mkdir -p virus
    cp ALL.fasta virus/ALL.fasta

    split_virus_db_files.py \\
        --yaml virus_input_config.yaml \\
        --fasta virus/ALL.fasta \\
        --outdir virus/ \\
        --description "${params.custom_db_virus_description}" \\
        --version "${params.custom_db_virus_version}" \\
        --blast-db blast_db \\
        --blast-prefix ALL \\
        --config-out virus.yaml \\
        --strict

    makeblastdb \
        -dbtype nucl \
        -in virus/ALL.fasta \
        -out virus/blast_db/ALL \
        -parse_seqids \
        -blastdb_version 5
    """
}


process custom_data_base_transfer {
    label "general"
    cpus 1
    publishDir (
        params.custom_db_outpath,
        mode: "copy",
        saveAs: {
            db
        }
    )
    input:
        path(db)
    output:
        path(db)
    """
    """
}


workflow custom_data_base {
    main:
        if(params.custom_db_contaminants_input_path) {
            contaminants_db = build_contaminants_db(params.custom_db_contaminants_input_path)
        } else {
            contaminants_db = Channel.empty()
        }

        if(params.custom_db_virus_yaml && params.custom_db_virus_fasta){
            virus_db = Channel.of([params.custom_db_virus_yaml, params.custom_db_virus_fasta])
            | build_virus_db
        } else {
            virus_db = Channel.empty()
        }

        centrifuge_db = Channel.empty()

        Channel.empty()
        | mix(
            contaminants_db,
            centrifuge_db,
            virus_db
        )
        | toList
        | flatMap
        | custom_data_base_transfer
}
