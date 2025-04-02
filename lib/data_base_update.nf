// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


nextflow.enable.dsl = 2


process db_update_get_config {
    label "general"
    cpus 1
    output:
        path("db.yaml")
    """
    wget -o download_config.log -O db.yaml ${params.update_db_config_url}
    """
}


def checksumDir(String directory) {
    return "find ${directory}/. -type f -exec sha256sum {} \\; | sort | awk '{print \$1}' | sha256sum | awk '{print \$1}'"
}


process check_download_necessary {
    label "general"
    cpus 1 
    input:
        tuple val(do_update), val(expected_checksum), path("db_path")
    output:
        path("dummy.txt"), optional: true
    """
    if [ "${do_update}" == "true" ]
    then
        if [[ -d "db_path" ]]
        then
            checksum_old_data=\$(${checksumDir("db_path")})
            if [ \$checksum_old_data != ${expected_checksum} ]
            then
                echo "update required because of different checksums (got \$checksum_old_data but expected ${expected_checksum})" > dummy.txt
            fi
        else
            echo "update required because the directory is missing" > dummy.txt
        fi
    fi
    """
}


process download_parts {
    maxParallel = 1
    label "download"
    cpus 1
    input:
        val(file_config)
    output:
        path("${file_config.name}")
    """
    wget -o dl_${file_config.name}.log -O ${file_config.name} ${file_config.url}

    checksum_new=\$(sha256sum ${file_config.name} | awk '{print \$1}')
    checksum_expected=${file_config.checksum}

    if [[ \$checksum_new != \$checksum_expected ]]
    then
        echo "Checksum mismatch!" >&2
        echo "Expected: \$checksum_expected" >&2
        echo "Computed: \$checksum_new" >&2
        exit 1
    fi
    """
}


process merge_parts_and_extract {
    label "general"
    cpus 1
    input:
        tuple val(db_config), path(parts), val(db_name)
    output:
        path("${db_name}")
    """
    # merge the files
    ls ${parts.join(' ')} | sort | xargs cat > merged.tar.xz

    # check the checksum
    checksum_new=\$(sha256sum merged.tar.xz | awk '{print \$1}')
    checksum_expected=${db_config.checksum_zipped}

    if [[ \$checksum_new != \$checksum_expected ]]
    then
        echo "Checksum mismatch for merged.tar.xz!" >&2
        echo "Expected: \$checksum_expected" >&2
        echo "Computed: \$checksum_new" >&2
        exit 1
    fi
    
    # extract
    tar -xf merged.tar.xz

    # check the directory checksum
    checksum_expected_directory=${db_config.checksum_directory}
    checksum_new_directory=\$(${checksumDir(db_name)})
    if [ "\$checksum_new_directory" != "\$checksum_expected_directory" ]
    then
        echo "Checksum mismatch for ${db_name}!" >&2
        echo "Expected: \$checksum_expected_directory" >&2
        echo "Computed: \$checksum_new_directory" >&2
        exit 1
    fi
    """
}


process data_base_transfer {
    label "general"
    cpus 1
    publishDir (
        params.database_defaults.base,
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
