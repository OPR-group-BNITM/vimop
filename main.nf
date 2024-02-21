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
    seqkit stats seqs.fastq > stats.txt
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
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}


// workflow module
workflow pipeline {
    take:
        samples
    main:
        // trimming
        trimmed = trim(samples.map{meta, reads, stats -> [meta, reads]})
        trim_stats = seqtk_stats(trimmed)
        trim_fastqc = fastqc(trimmed)

        to_classify = trimmed.combine(Channel.of(params.virus_db).combine(Channel.of('viral', 'hpvc')))
        classification = classify_centrifuge(to_classify)

        // define output
        ch_to_publish = Channel.empty()
        | mix(
            trim_stats | map {meta, stats -> [stats, "$meta.alias/trim_stats"]},
            trim_fastqc | map {meta, fastqc -> [fastqc, "$meta.alias/trim_stats"]},
            classification | map {meta, classification, report, kraken, html -> [classification, "$meta.alias/classification"]},
            classification | map {meta, classification, report, kraken, html -> [report, "$meta.alias/classification"]},
            classification | map {meta, classification, report, kraken, html -> [kraken, "$meta.alias/classification"]},
            classification | map {meta, classification, report, kraken, html -> [html, "$meta.alias/classification"]},
        )

        // TODO cleaning

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
