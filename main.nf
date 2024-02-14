// workflow module
    workflow pipeline {
        take:
            samples
            reference
            blastdb
            nextclade_data
        main:
            alignment = alignReads(samples.map{ meta, reads, stats -> [ meta,reads ] }, reference)
            coverage = coverageCalc(alignment.alignments)

            // do crude downsampling
            if (params.rbk){
                println("RBK data - NOT Downsampling!!!")
                downsample = alignment
            } else if (params.downsample != null){
                println("Downsampling!!!")
                downsample = downSample(alignment.alignments, reference)
            } else {
                println("NOT Downsampling!!!")
                downsample = alignment
            }

            if(params.medaka_consensus_model) {
                log.warn "Overriding Medaka consensus model with ${params.medaka_consensus_model}."
                medaka_consensus_model = Channel.fromList([params.medaka_consensus_model])
            }
            else {
                lookup_table = Channel.fromPath("${projectDir}/data/medaka_models.tsv", checkIfExists: true)
                medaka_consensus_model = lookup_medaka_consensus_model(lookup_table, params.basecaller_cfg)
            }

            bams_for_calling = downsample.alignments.combine(medaka_consensus_model)

            variants = medakaVariants(bams_for_calling, reference)

            for_draft = variants.join(coverage)

            draft = makeConsensus(for_draft, reference)
            flu_type = typeFlu(draft, blastdb)

            processed_type = processType(flu_type.typing)


            nextclade_prep = prepNextclade(
                processed_type.join(draft, remainder: true),
                nextclade_data
            )

            nextclade_datasets = nextclade_prep
            | map { file(it.resolve("**"), type: "file") }
            | flatten
            | map { [it.parent.name, it] }
            | groupTuple


            nextclade_result = nextclade(nextclade_datasets)

            software_versions = getVersions()
            software_versions = software_versions.mix(flu_type.version.first()) | collectFile()
            workflow_params = getParams()

            // output_alignments = alignment.alignments.map{ it -> return tuple(it[1], it[2]) }

            // get all the per sample results together
            ch_per_sample_results = samples
            | join(coverage)
            | join(flu_type.typing) 
            | join(processed_type)

            // collect results into a directory for the sample directory to avoid collisions
            ch_results_for_report = ch_per_sample_results
            | map {
                meta = it[0]
                rest = it[1..-1]
                [meta, rest, meta.alias]
            }
            | collectFilesInDir
            | map { meta, dirname -> dirname }



            report = makeReport(
                samples.map{it -> it[0]}.toList(),
                samples.map{it -> it[2].resolve("per-read-stats.tsv.gz")}.toList(),
                ch_results_for_report | collect,
                nextclade_result.map{it -> it[1]}.ifEmpty(OPTIONAL_FILE).collect(),
                nextclade_data,
                software_versions.collect(),
                workflow_params
            )

            // create channel with files to publish; the channel will have the shape `[file,
            // name of sub-dir to be published in]`.

            ch_to_publish = Channel.empty()
            | mix(
                software_versions | map { [it, null] },
                workflow_params | map { [it, null] },
                report | map { [it[0] , null] },
                report | map { [it[1], null] },
                alignment.alignments
                | map { meta, bam, bai -> [bam, "$meta.alias/alignments"] },
                variants
                | map { meta, vcf -> [vcf, "$meta.alias/variants"]},
                coverage
                | map { meta, depth -> [depth, "$meta.alias/coverage"]},
                draft
                | map { meta, fa -> [fa, "$meta.alias/consensus"]},
                flu_type.typing
                | map { meta, json -> [json, "$meta.alias/typing"]}
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
        "sample_sheet": params.sample_sheet]
    )


    // TODO

//get reference
    if (params.reference == null){
        params.remove('reference')
        params._reference = projectDir.resolve("./data/primer_schemes/V1/consensus_irma.fasta").toString()
    } else {
        params._reference = file(params.reference, type: "file", checkIfExists:true).toString()
        params.remove('reference')
    }

    //get reference
    if (params.blastdb == null){
        params.remove('blastdb')
        params._blastdb = projectDir.resolve("./data/primer_schemes/V1/blastdb").toString()
    } else {
        params._blastdb = file(params.reference, type: "file", checkIfExists:true).toString()
        params.remove('blastdb')
    }

    nextclade_data = projectDir.resolve("./data/nextclade.csv").toString()

    pipeline(samples, params._reference, params._blastdb, nextclade_data)
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
