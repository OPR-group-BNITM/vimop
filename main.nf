import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'


process trim {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("demultiplexed.fastq.gz")
    output:
        tuple val(meta), path("trimmed.fastq")
    """
    seqtk trimfq -b 30 -e 30 demultiplexed.fastq.gz > trimmed.fastq
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
    contaminants=(${contaminants.join(" ")})
    fn_input=seqs.fastq

    echo "step\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > stats_tmp.tsv
    echo -n "nofilter\t" >> stats_tmp.tsv
    seqkit stats -T \$fn_input | tail -n 1 | awk '{print \$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8}' >> stats_tmp.tsv

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

        echo -n "\${c}\t" >> stats_tmp.tsv
        seqkit stats -T \$fn_out | tail -n 1 | awk '{print \$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8}' >> stats_tmp.tsv

        fn_input=\$fn_out
    done
    cp \$fn_input filtered.fastq
    mv stats_tmp.tsv stats.tsv
    """
}


process filter_virus_target {
    label "opr_target_mapping"
    cpus 4
    input:
        tuple val(meta), path('seqs.fastq'), path(db_path)
    output:
        tuple val(meta), path('filtered.fastq')
    """
    if [[ -f ${db_path}/${meta.mapping_target}.mmi ]]; then
        fname_ref_index=${db_path}/${meta.mapping_target}.mmi
    elif [[ -f ${db_path}/${meta.mapping_target}.fasta ]]; then
        fname_ref_index=${db_path}/${meta.mapping_target}.fasta
    else
        echo "No file ${db_path}/${meta.mapping_target}.fasta or ${db_path}/${meta.mapping_target}.mmi found as reference file to map" >&2
        exit 1
    fi

    minimap2 \
        -x map-ont \
        -a \$fname_ref_index \
        -t ${task.cpus} \
        seqs.fastq \
    | samtools fastq \
        -F 4 \
        --reference reference.fasta > filtered.fastq
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
    subsampling_sizes=(${subsampling_sizes.join(" ")})
    min_readlengths=(${min_readlengths.join(" ")})
    min_overlaps=(${min_overlaps.join(" ")})

    echo "step\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > stats.tsv
    echo -n "all_${meta.mapping_target}\t" >> stats.tsv
    seqkit stats -T seqs.fastq | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv

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

        set +e
        canu \
            -nanopore-raw \$fname_subsampled \
            -fast \
            -d \$outdir \
            -p \$prefix \
            genomeSize=$genome_size \
            minReadLength=\$min_readlen \
            minOverlapLength=\$min_overlap \
            corOutCoverage=$cor_out_coverage
        canu_status=\$?
        set -e
        if [[ \$canu_status -eq 0 ]] && [[ -f \$fname_contigs ]] && [[ -s \$fname_contigs ]]; then
            mv \$fname_contigs assemblies.fasta
            break
        else
            continue
        fi
    done

    # creates empty file if not exist
    touch assemblies.fasta
    """
}


process prepare_blast_search {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("contigs.fasta")
    output:
        tuple val(meta), path("sorted-contigs.fasta")
    """
    seqkit sort --by-length --reverse contigs.fasta \
    | seqkit rename > sorted-contigs.fasta
    """
}


process contigs_stats {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("sorted-contigs.fasta")
    output:
        tuple val(meta), path("stats.tsv")
    """
    seqkit stats sorted-contigs.fasta > stats.tsv
    """
}


process contigs_readlengths {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("sorted-contigs.fasta")
    output:
        tuple val(meta), path("readlengths.txt")
    """
    readlength.sh qin=33 in=sorted-contigs.fasta bin=1 out=readlengths.txt
    """
}


process blast {
    label "opr_general"
    cpus 4
    input:
        tuple val(meta), path("sorted-contigs.fasta"), path(db_path), val(target)
    output:
        tuple val(meta), path("blast-results.xml")
    """
    blastn \
        -num_threads ${task.cpus} \
        -db ${db_path}/${target} \
        -query sorted-contigs.fasta \
        -out blast-results.xml \
        -outfmt 5
    """
}


process extract_blasthits {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("blast-results.xml")
    output:
        tuple val(meta), path("blast-hits.csv")
    """
    #!/usr/bin/env python
    import pandas as pd
    from xml.etree import ElementTree
    rows = []
    try:
        root = ElementTree.parse('blast-results.xml').getroot()
        for hit in root.findall(".//Hit[Hit_num='1']"):
            try:
                accession = hit.find('Hit_accession').text
                description = hit.find('Hit_def').text.replace(',', '').replace(';', '')
                length = int(hit.find('Hit_len').text)
                bit_score = float(hit.find(".//Hsp_bit-score").text)
                rows.append(['${meta.alias}', '${meta.mapping_target}', accession, description, length, bit_score])
            except (ValueError, AttributeError):
                continue
    except ElementTree.ParseError:
        # empty file with no blast hits
        pass
    hits = pd.DataFrame(rows, columns=['sample', 'mapping_target', 'ref', 'def', 'length', 'bitscore'])
    hits_bitscore_summed = hits.groupby(['sample', 'mapping_target', 'ref','def','length']).aggregate({'bitscore': 'sum'}).reset_index()
    hits_bitscore_summed.to_csv('blast-hits.csv', header=True, index=False)
    """
}


process get_ref_fasta {
    label "opr_general"
    cpus 1
    input:
        tuple val(ref_id), path(db_path), val(db_name)
    output:
        tuple val(ref_id), path("ref.fasta")
    """
    blastdbcmd -entry ${ref_id} -db ${db_path}/${db_name} -out ref.fasta
    """
}


process map_to_ref {
    label "opr_map"
    cpus 8
    input:
        tuple val(meta), path("trimmed.fastq"), path("ref.fasta")
    output:
        tuple val(meta), path("ref.fasta"), path("sorted.bam"), path("sorted.bam.bai")
    """
    minimap2 -ax map-ont ref.fasta trimmed.fastq -t ${task.cpus} | samtools view -b -o mapped.bam
    samtools sort --threads ${task.cpus} -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}


process medaka_consensus {
    label "opr_medaka"
    cpus 2
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            path("model.hdf5"),
            val(min_depth)
    output:
        tuple val(meta), path('consensus.fasta')
    """
    medaka consensus sorted.bam consensus.hdf --threads ${task.cpus} --model model.hdf5
    medaka stitch consensus.hdf ref.fasta consensus.fasta --fill_char N --min_depth ${min_depth}
    """
    // todo use --qualities in stitch. then filter mask based on the quality? and export the consensus as fastq or fasta?
    // seqtk seq -q20 -n N input.fastq > output.fastq
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
        trimmed = samples | map{ meta, reads, stats -> [meta, reads] } | trim
        trim_stats = trimmed | seqtk_stats
        trim_fastqc = trimmed | fastqc

        // metagenomic read classification with centrifuge
        classification = trimmed
        | combine(Channel.of(params.virus_db))
        | combine(Channel.from(params.centrifuge_classification_libraries))
        | classify_centrifuge

        // contaminant filtering
        db_files = params.contamination_filters.collect{ filter -> params.contamination_databases[filter] }
        cleaned = trimmed
        | map{ meta, reads -> [meta, reads, db_files, params.contamination_filters] }
        | filter_contaminants

        // mapping reads to given virus targets to filter them
        mapped_to_virus_target = cleaned
        | combine(Channel.from(params.targets))
        | map{ meta, reads, stats, target -> [meta + ["mapping_target": target], reads, params.virus_db]}
        | filter_virus_target

        // assembly to get queries to find references
        def assembly_params = [
            params.assembly_parameters.collect { it.n_reads },
            params.assembly_parameters.collect { it.min_readlen },
            params.assembly_parameters.collect { it.min_overlap },
            params.assembly_target_genome_size,
            params.assembly_cor_out_coverage
        ]

        to_assemble_targeted = mapped_to_virus_target
        | map { meta, reads -> [meta, reads] + assembly_params }

        if (params.assemble_notarget) {
            to_assemble_notarget = cleaned
            | map { meta, reads, stats -> [meta + ["mapping_target": "no-target"], reads] + assembly_params }
        } else {
            to_assemble_notarget = Channel.empty()
        }
        assemblies = to_assemble_notarget | mix(to_assemble_targeted) | assemble_canu

        // TODO: if no draft assembly -> longest reads? (or are the longest reads automatically in the assemblies?)

        // database search for references using blast
        blast_queries = assemblies | map { meta, contigs, stats -> [meta, contigs] } | prepare_blast_search
        blast_query_stats = blast_queries | contigs_stats
        blast_query_lengths = blast_queries | contigs_readlengths
        blast_hits = blast_queries
        | map { meta, contigs -> [meta, contigs, params.virus_db, params.all_target] }
        | blast
        | extract_blasthits

        // get unique mapping targets from blast hits stored in .csv files
        sample_ref = blast_hits
        | flatMap {meta, hits -> hits.readLines().drop(1).collect { line -> tuple(meta.alias, line.split(',')[2]) }}
        | unique

        // TODO:
        // add manual references (sample + ref_accession) !!!!!!!

        // get the targets
        ref_seqs = sample_ref
        | map { samplename, ref_id -> ref_id }
        | unique
        | map { ref_id -> [ref_id, params.virus_db, params.all_target]}
        | get_ref_fasta

        reads_and_ref = sample_ref
        | combine(ref_seqs)
        | filter { samplename, ref_id_1, ref_id_2, ref_seq -> ref_id_1 == ref_id_2 }
        | map { samplename, ref_id_1, ref_id_2, ref_seq -> [samplename, ref_id_1, ref_seq] }
        | combine(trimmed)
        | filter { samplename, ref_id, ref_seq, meta, reads -> samplename == meta.alias }
        | map { samplename, ref_id, ref_seq, meta, reads -> [meta + ["consensus_target": ref_id], reads, ref_seq] }

        // Map against references
        mapped_to_ref = reads_and_ref | map_to_ref

        // build the consensus sequences
        consensi = mapped_to_ref
        | map { meta, ref, bam, bai -> [meta, ref, bam, bai, params.medaka_consensus_model, params.consensus_min_depth] }
        | medaka_consensus

        // TODO: incorporate minimum coverage?
        // TODO: create vcf files
        // TODO: export more consensus files?
        // TODO: export vcf files

        // TODO: assemble report - .csv (pandas)
        // TODO: create html reports

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
            assemblies | map { meta, assemblies, stats -> [stats, "$meta.alias/assembly/$meta.mapping_target", null] },
            blast_hits | map { meta, hits -> [hits, "$meta.alias/blast-hits/$meta.mapping_target", null] },
            blast_query_stats | map { meta, stats -> [stats, "$meta.alias/blast-hits/$meta.mapping_target", null] },
            blast_query_lengths | map { meta, lengths -> [lengths, "$meta.alias/blast-hits/$meta.mapping_target", null] },
            mapped_to_ref | map { meta, ref, bam, bai -> [ref, "$meta.alias/consensus/mapping", "${meta.consensus_target}.fasta"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bam, "$meta.alias/consensus/mapping", "${meta.consensus_target}.bam"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bai, "$meta.alias/consensus/mapping", "${meta.consensus_target}.bam.bai"] },
            consensi | map { meta, consensus -> [consensus, "$meta.alias/consensus/consensus/${meta.consensus_target}", null] },
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
