import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'


process trim {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("demultiplexed.fastq")
    output:
        tuple val(meta), path("trimmed.fastq")
    """
    seqtk trimfq -b 30 -e 30 demultiplexed.fastq > trimmed.fastq
    """
}


process seqtk_stats {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("seqs.fastq")
    output:
        tuple val(meta), path("stats.txt")
    """
    seqkit stats -T seqs.fastq > stats.txt
    """
}


process fastqc {
    label "opr_general"
    cpus 4
    input:
        tuple val(meta), path("seqs.fastq")
    output:
        tuple val(meta), path("fastqc")
    """
    mkdir -p fastqc
    fastqc -q -t ${task.cpus} -o fastqc seqs.fastq 
    """
}


process classify_centrifuge {
    label "opr_centrifuge"
    cpus 1
    input:
        tuple val(meta), path("seqs.fastq"), path(db_path), val(target_db)
    output:
        tuple val(meta),
              path("classification-${target_db}.tsv"),
              path("classification-report-${target_db}.tsv"),
              path("classification-kraken-${target_db}.tsv"),
              path("classification-${target_db}.html")
    """
    centrifuge \
        --mm \
        -x ${db_path}/${target_db} \
        -U seqs.fastq \
        --report-file classification-report-${target_db}.tsv \
        -S classification-${target_db}.tsv
    centrifuge-kreport \
        -x ${db_path}/${target_db} classification-${target_db}.tsv > classification-kraken-${target_db}.tsv
    ktImportTaxonomy \
        -tax ${db_path}/taxonomy \
        -m 3 -t 5 \
        classification-kraken-${target_db}.tsv \
        -o classification-${target_db}.html
    """
}


process filter_contaminants {
    label "opr_clean"
    cpus 4
    input:
        tuple val(meta), path('seqs.fastq'), path('db_*.fna.gz'), val(contaminants)
    output:
        tuple val(meta), path('filtered.fastq'), path('stats.tsv')
    """
    #!/usr/bin/env bash

    contaminants=(${contaminants.join(" ")})
    fn_input=seqs.fastq

    echo "step\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > stats.tsv
    echo -n "nofilter\t" >> stats.tsv
    seqkit stats -T \$fn_input | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv

    for i in "\${!contaminants[@]}"
    do
        c=\${contaminants[\$i]}
        db=db_\$((i+1)).fna.gz
        fn_out=filtered_no_\${c}.fastq

        fn_sam=filter_\${c}.sam

        minimap2 \
            -x map-ont \
            -a \$db \
            -t ${task.cpus} \
            \$fn_input \
        | samtools fastq \
            -f 4 \
            --reference \$db \
            -1 \${fn_out} -2 \${fn_out} -0 \${fn_out} -s \${fn_out} -n

        echo -n "\${c}\t" >> stats.tsv
        seqkit stats -T \$fn_out | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv

        fn_input=\$fn_out
    done
    mv \$fn_input filtered.fastq
    """
}


process filter_virus_target {
    label "opr_target_mapping"
    cpus 6
    input:
        tuple val(meta), path('seqs.fastq'), path('reference.mmi'), path('reference.fasta')
    output:
        tuple val(meta), path('filtered.fastq'), path('stats.tsv')
    """
    #!/usr/bin/env bash

    minimap2 \
        -x map-ont \
        -a reference.mmi \
        -t 4 \
        seqs.fastq \
    | samtools fastq \
        --threads 2 \
        -F 4 \
        --reference reference.fasta > filtered.fastq
    
    seqkit stats -T filtered.fastq > stats.tsv
    """
}


process assemble_canu {
    label "opr_draft_assembly"
    cpus 8
    input:
        tuple val(meta),
            path("seqs.fastq"),
            val(subsampling_sizes),
            val(min_readlengths),
            val(min_overlaps),
            val(genome_size),
            val(cor_out_coverage)
    output:
        tuple val(meta), path("assemblies.fasta"), path("stats.tsv")
    """
    #!/usr/bin/env bash

    subsampling_sizes=(${subsampling_sizes.join(" ")})
    min_readlengths=(${min_readlengths.join(" ")})
    min_overlaps=(${min_overlaps.join(" ")})

    echo "step\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > stats.tsv

    for i in "\${!subsampling_sizes[@]}"
    do
        n_reads=\${subsampling_sizes[\$i]}
        min_readlen=\${min_readlengths[\$i]}
        min_overlap=\${min_overlaps[\$i]}

        # check if the subsampled file exists, else subsample
        fname_subsampled=seqs_\${n_reads}.fastq
        if [ ! -f \$fname_subsampled ]; then
            reformat.sh qin=33 ow=t samplereadstarget=\${n_reads} in=seqs.fastq out=\${fname_subsampled}
            echo -n "subsample_\${n_reads}\t" >> stats.tsv
            seqkit stats -T \$fname_subsampled | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv
        fi

        outdir=attempt_\$i
        prefix=attempt_\$i
        fname_contigs=\${outdir}/\${prefix}.contigs.fasta

        canu \
            -nanopore-raw \$fname_subsampled \
            -fast \
            -d \$outdir \
            -p \$prefix \
            genomeSize=$genome_size \
            minReadLength=min_readlen \
            minOverlapLength=min_overlap \
            corOutCoverage=$cor_out_coverage

        canu_status=\$?

        if [[ \$canu_status -eq 0 ]] && [[ -f \$fname_contigs ]] && [[ -s \$fname_contigs ]]; then
            mv \$fname_contigs assemblies.fasta
            break
        fi
    done

    # creates empty file if not exist
    touch assemblies.fasta
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wf_common"
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { 
            dirname ? (fname_out ? "$dirname/$fname_out" : "$dirname/$fname_in") : fname_in
        }
    )
    input:
        tuple path(fname_in), val(dirname), val(fname_out)
    output:
        path(fname_in)
    """
    """
}


// workflow module
workflow pipeline {
    take:
        samples
    main:
        // trimming
        trimmed = trim(samples.map{ meta, reads, stats -> [meta, reads] })
        trim_stats = seqtk_stats(trimmed)
        trim_fastqc = fastqc(trimmed)

        // TODO: add classification databases to configs!
        // metagenomic read classification with centrifuge
        // to_classify = trimmed.combine(Channel.of(params.virus_db).combine(Channel.of('viral', 'hpvc')))
        to_classify = trimmed.combine(Channel.of(params.virus_db).combine(Channel.of('viral')))
        classification = classify_centrifuge(to_classify)

        // contaminant filtering
        db_files = params.contamination_filters.collect{ filter -> params.contamination_databases[filter] }
        cleaned = filter_contaminants(
            trimmed.map{ meta, reads -> [meta, reads, db_files, params.contamination_filters] }
        )

        // TODO: create the mapping indices? - precompute for the whole DB!
        // - get the targets from the params, also get the minimap version from the parameters?
        // TODO: mapping

        def assembly_params = [
            params.assembly_parameters.collect { it.n_reads },
            params.assembly_parameters.collect { it.min_readlen },
            params.assembly_parameters.collect { it.min_overlap },
            params.assembly_target_genome_size,
            params.assembly_cor_out_coverage
        ]

        // TODO: only if assemble without target
        // TODO: add target-mapped results
        to_assemble_targeted = Channel.empty()

        if (params.assemble_notarget) {
            to_assemble_notarget = cleaned.map { meta, reads, stats -> [meta + ["mapping_target": "no-target"], reads] + assembly_params }
        } else {
            to_assemble_notarget = Channel.empty()
        }
        to_assemble = to_assemble_notarget.mix(to_assemble_targeted)

        assemblies = assemble_canu(to_assemble)

        // TODO: get blast candidates

        // TODO: blast

        // TODO: get references from blast

        // TODO: map reads against references

        // TODO: create consensus

        // TODO: assemble report

        // define output
        ch_to_publish = Channel.empty()
        | mix(
            trim_stats | map { meta, stats -> [stats, "$meta.alias/trim_stats", null] },
            trim_fastqc | map { meta, fastqc -> [fastqc, "$meta.alias/trim_stats", null] },
            classification | map { meta, classification, report, kraken, html -> [classification, "$meta.alias/classification", null] },
            classification | map { meta, classification, report, kraken, html -> [report, "$meta.alias/classification", null] },
            classification | map { meta, classification, report, kraken, html -> [kraken, "$meta.alias/classification", null] },
            classification | map { meta, classification, report, kraken, html -> [html, "$meta.alias/classification", null] },
            cleaned | map { meta, reads, stats -> [stats, "$meta.alias/clean", null] },
            assemblies | map { meta, assemblies, stats -> [assemblies, "$meta.alias/assembly/$meta.mapping_target", null] },
            assemblies | map { meta, assemblies, stats -> [stats, "$meta.alias/assembly/$meta.mapping_target", null] }
        )

    emit:
        results = ch_to_publish
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)
    samples = fastq_ingress([
        "input": params.fastq,
        "stats": true,
        "sample_sheet": params.sample_sheet
    ])
    // samples.view {
    //     sample ->
    //     println "${sample[0].barcode} ${sample[0].type} ${sample[0].run_ids} ${sample[0].alias}\n${sample[1]}\n${sample[2]}"
    // }
    pipeline(samples)
    pipeline.out.results
    | toList
    | flatMap
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
