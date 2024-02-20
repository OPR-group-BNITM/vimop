import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'


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
