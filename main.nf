import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'


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
        reference
        blastdb
        nextclade_data
    main:
        // TODO
        // alignment = alignReads(samples.map{ meta, reads, stats -> [ meta,reads ] }, reference)

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

    // print the output samples
    samples.view {
        sample->
        println "${sample[0].barcode} ${sample[0].type} ${sample[0].run_ids} ${sample[0].alias}\n${sample[1]}\n${sample[2]}"
    }

    // TODO
    // Run a workflow
    // publish data to output directory
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
