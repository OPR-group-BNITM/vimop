import groovy.json.JsonBuilder
import DatabaseInput
import SystemRequirements

nextflow.enable.dsl = 2


include { fastq_ingress } from './lib/ingress'
include {
    lengths_and_qualities as lengths_and_qualities_trimmed;
    lengths_and_qualities as lengths_and_qualities_cleaned;
    trim;
    classify_centrifuge;
    filter_contaminants;
    filter_virus_target;
    assemble_canu;
    prepare_blast_search;
    canu_contig_info;
    blast;
    extract_blasthits;
    get_ref_fasta;
    map_to_ref;
    calc_coverage;
    medaka_consensus;
    simple_consensus;
    compute_mapping_stats;
    concat_mapping_stats;
    collect_reference_info;
    sample_report;
    get_best_consensus_files;
    output;
} from './lib/processes'


// workflow module
workflow pipeline {
    take:
        samples
    main:
        db_config = new DatabaseInput(params)

        // trimming
        trimmed = samples
        | map{ meta, reads, stats -> [meta, reads] }
        | trim

        // metagenomic read classification with centrifuge
        classification = trimmed
        | combine(Channel.of(db_config.classificationDir))
        | combine(Channel.from(db_config.classificationLibraries))
        | classify_centrifuge

        // contaminant filtering
        cleaned = trimmed
        | map{ meta, reads -> [meta, reads, db_config.contaminationFilterFiles, db_config.contaminationFilters] }
        | filter_contaminants

        // get readstats
        lenquals_trim = trimmed
        | map {meta, reads -> [meta.alias, reads]}
        | lengths_and_qualities_trimmed
    
        lenquals_clean = cleaned.reads
        | map {meta, reads -> [meta.alias, reads]}
        | lengths_and_qualities_cleaned

        // mapping reads to given virus targets to filter them
        mapped_to_virus_target = cleaned.reads
        | combine(Channel.from(db_config.virusTargets))
        | map{ meta, reads, target -> [meta + ["mapping_target": target.target], reads, target.path] }
        | filter_virus_target

        // assembly to get queries to find references
        def assembly_params = [
            params.canu_min_read_length,
            params.canu_min_overlap_length,
            params.canu_read_sampling_bias,
            params.canu_genome_size,
            params.canu_cor_out_coverage,
            params.canu_stop_on_low_coverage,
            params.canu_min_input_coverage
        ]

        to_assemble_targeted = mapped_to_virus_target
        | map { meta, reads -> [meta, reads] + assembly_params }

        if (params.assemble_notarget) {
            to_assemble_notarget = cleaned.reads
            | map { meta, reads -> [meta + ["mapping_target": "no-target"], reads] + assembly_params }
        } else {
            to_assemble_notarget = Channel.empty()
        }
        assemblies = to_assemble_notarget
        | mix(to_assemble_targeted)
        | assemble_canu

        // database search for references using blast
        blast_queries = assemblies.contigs
        | pop_bubbles
        | prepare_blast_search

        collected_contigs_infos = blast_queries
        | canu_contig_info
        | map { meta, contig_info -> [meta.alias, contig_info]}
        | groupTuple(by: 0)

        blast_hits = blast_queries
        | map { meta, contigs -> [meta, contigs, db_config.blastDir, db_config.blastPrefix] }
        | blast
        | extract_blasthits

        // get mapping targets from blast hits stored in .csv files
        sample_ref = blast_hits
        | flatMap {meta, hits -> hits.readLines().drop(1).collect { line -> tuple(meta.alias, meta.mapping_target, line.split(',')[1]) }}
        | map {samplename, target_filter, blast_hit -> [samplename, blast_hit]}
        | unique

        // get custom accessions for references given by user
        custom_refs = Channel.from(params.custom_sample_ref)
        // add custom references to the sample_ref channel
        sample_ref = sample_ref
        | concat(custom_refs)
        | unique

        // get the targets
        ref_seqs = sample_ref
        | map { samplename, ref_id -> ref_id }
        | unique
        | map { ref_id -> [ref_id, db_config.blastDir, db_config.blastPrefix] }
        | get_ref_fasta

        // add custom reference files here
        custom_sample_ref_files = Channel.from(params.custom_sample_ref_files)

        // read custom reference files
        custom_sample_ref_id_seqs = custom_sample_ref_files
        | map { sample, file_path ->
            lines = file(file_path).readLines()
            matches = lines.findAll { line -> line =~ /^>/ } // get Fasta headers
                        .collectMany { line -> (line =~ /(?<=>)[^\s]+/).findAll() } // Extract entry_ids from the header
            // Assert that there is only one sequence per file
            assert matches.size() == 1 : "ERROR: Custom reference file must contain only one sequence: ${file_path}"
            custom_sample_ref_ids = matches.collect { match -> [sample, match] }.flatten() // Pair each match with the sample
            custom_sample_ref_seqs = matches.collect { match -> [match, file_path] }.flatten() // Pair each match with the file_path
            [custom_sample_ref_ids, custom_sample_ref_seqs]
        }

        custom_ref_seqs = custom_sample_ref_id_seqs // [ref_id, path_to_ref_seq]
        | flatMap { custom_sample_ref_ids, custom_sample_ref_seqs -> [custom_sample_ref_seqs] }

        custom_sample_ref_ids = custom_sample_ref_id_seqs // [sample, ref_id]
        | flatMap { custom_sample_ref_ids, custom_sample_ref_seqs -> [custom_sample_ref_ids] }

        // extend the reference sequences with the custom ones
        extended_ref_seqs = ref_seqs
        | concat(custom_ref_seqs)
        | unique

        extended_sample_ref = sample_ref
        | concat(custom_sample_ref_ids)
        | unique

        reads_and_ref = extended_sample_ref
        | combine(extended_ref_seqs)
        | filter { samplename, ref_id_1, ref_id_2, ref_seq -> ref_id_1 == ref_id_2 }
        | map { samplename, ref_id_1, ref_id_2, ref_seq -> [samplename, ref_id_1, ref_seq] }
        | combine(trimmed)
        | filter { samplename, ref_id, ref_seq, meta, reads -> samplename == meta.alias }
        | map { samplename, ref_id, ref_seq, meta, reads -> [meta + ["consensus_target": ref_id], reads, ref_seq] }

        // Map against references
        mapped_to_ref = reads_and_ref
        | map_to_ref

        // get coverages
        coverage = mapped_to_ref
        | map { meta, ref, bam, bai -> [meta, bam, bai]}
        | calc_coverage

        if (params.consensus_method == 'simple') {
            consensi = mapped_to_ref
            | map { meta, ref, bam, bai -> [meta, ref, bam, bai, params.consensus_min_depth, params.consensus_min_share] }
            | simple_consensus
        } else if (params.consensus_method == 'medaka') {
            consensi = mapped_to_ref
            | map { meta, ref, bam, bai -> [meta, ref, bam, bai, params.medaka_consensus_model, params.consensus_min_depth] }
            | medaka_consensus
        }

        reference_info = extended_ref_seqs
        | map {refid, refseq -> refseq}
        | collect
        | collect_reference_info

        // TODO: create & export vcf files?

        // compute the stats for the reads mapped to the different targets and build a big table.
        collected_mapping_stats = mapped_to_ref
        | map {meta, ref, bam, bai -> [meta.alias, meta.consensus_target, meta, ref, bam, bai]}
        | join(consensi | map {meta, cons -> [meta.alias, meta.consensus_target, cons]}, by: [0, 1])
        | join(coverage | map {meta, cov -> [meta.alias, meta.consensus_target, cov]}, by: [0, 1])
        | map {samplename, target_name, meta, ref, bam, bai, cons, cov -> [meta, ref, bam, bai, cons, cov]}
        | compute_mapping_stats
        | map {meta, stats -> [meta.alias, stats]}
        | groupTuple(by: 0)
        | concat_mapping_stats

        // Create the report
        assembly_modes = params.assemble_notarget ? params.targets + ['no-target'] : params.targets

        collected_assembly_stats = assemblies.stats
        | map {meta, stats -> [meta.alias, stats]}
        | groupTuple(by: 0)

        collected_blast_hits = blast_hits
        | map {meta, hits -> [meta.alias, hits]}
        | groupTuple(by: 0)

        sample_results = cleaned.stats
        | map {samplename, clean_stats -> [samplename, clean_stats, assembly_modes]}
        | join(collected_assembly_stats, by: 0)
        | join(collected_blast_hits, by: 0)
        | join(collected_mapping_stats, by: 0)
        | join(collected_contigs_infos, by: 0)
        | join(lenquals_trim, by: 0)
        | join(lenquals_clean, by: 0)
        | combine(Channel.of(db_config.virusConfigFileName))
        | combine(reference_info)
        | sample_report

        collected_consensi = consensi
        | map {meta, consensus -> [meta.alias, consensus]}
        | groupTuple(by: 0)

        best_consensi = sample_results.consensus_stats
        | join(collected_consensi, by: 0)
        | get_best_consensus_files

        // define output
        ch_to_publish = Channel.empty()
        | mix(
            // centrifuge classification
            classification | map { meta, classification, report, kraken, html -> [classification, "$meta.alias/classification", null] },
            classification | map { meta, classification, report, kraken, html -> [report, "$meta.alias/classification", null] },
            classification | map { meta, classification, report, kraken, html -> [kraken, "$meta.alias/classification", null] },
            classification | map { meta, classification, report, kraken, html -> [html, "$meta.alias/classification", null] },
            // contigs
            assemblies.contigs | map { meta, contigs -> [contigs, "$meta.alias/assembly", "${meta.mapping_target}.contigs.fasta"] },
            // consensus
            mapped_to_ref | map { meta, ref, bam, bai -> [ref, "$meta.alias/consensus", "${meta.consensus_target}.reference.fasta"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bam, "$meta.alias/consensus", "${meta.consensus_target}.reads.bam"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bai, "$meta.alias/consensus", "${meta.consensus_target}.reads.bam.bai"] },
            consensi | map { meta, consensus -> [consensus, "$meta.alias/consensus", "${meta.consensus_target}.consensus.fasta"] },
            coverage | map { meta, coverage -> [coverage, "$meta.alias/consensus", "${meta.consensus_target}.depth.txt"] },
            // selected consensi
            best_consensi | map { alias, consensus_dir -> [consensus_dir, "$alias", "selected_consensus"] },
            // report and tables
            sample_results.report | map { alias, report -> [report, "$alias", "report_${alias}.html"] },
            sample_results.read_stats | map { alias, read_stats -> [read_stats, "$alias/tables", "reads.tsv"] },
            sample_results.contig_stats | map { alias, contig_stats -> [contig_stats, "$alias/tables", "contigs.tsv"] },
            sample_results.consensus_stats | map { alias, consensus_stats -> [consensus_stats, "$alias/tables", "consensus.tsv"] },
            // advanced output
            cleaned.reads | filter { params.output_cleaned_reads } | map { meta, reads -> [reads, "$meta.alias/clean", null] }
        )

    emit:
        results = ch_to_publish
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    def workDir = session.workDir.toString()
    new SystemRequirements(true).checkSystemRequirements(
        params.min_disk_space_work_gb,
        params.min_disk_space_out_gb,
        params.min_ram_gb,
        params.min_cpus,
        params.out_dir,
        workDir
    )

    samples = fastq_ingress([
        "input": params.fastq,
        "stats": true
    ])
    pipeline(samples)
    pipeline.out.results
    | toList
    | flatMap
    | output
}


workflow.onComplete {
    File outputFile = new File("${params.out_dir}/params.json")
    def json = new JsonBuilder(params)
    outputFile.withWriter('UTF-8') { writer -> writer.write(json.toPrettyString()) }
    // notify that the workflow is complete
    Pinguscript.ping_complete(nextflow, workflow, params)
}


workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
