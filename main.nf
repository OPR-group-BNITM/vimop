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
    seqtk trimfq -b ${params.trim_length} -e ${params.trim_length} demultiplexed.fastq.gz > trimmed.fastq
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
        tuple val(meta), path('filtered.fastq'), path("clean_stats.tsv")
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
    mv stats_tmp.tsv clean_stats.tsv
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
    memory '30.GB'
    input:
        tuple val(meta),
            path("seqs.fastq"),
            val(subsampling_sizes),
            val(min_readlengths),
            val(min_overlaps),
            val(genome_size),
            val(cor_out_coverage)
    output:
        tuple val(meta),
            path("assemblies.fasta"),
            path("assembly_stats_${meta.mapping_target}.tsv")
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
        fname_minlen=minlen_\${min_readlen}.fastq
        if [ ! -f \$fname_minlen ]; then
            reformat.sh qin=33 ow=t in=seqs.fastq out=\$fname_minlen minlen=\$min_readlen
            echo -n "minlen_\${min_readlen}\t" >> stats.tsv
            seqkit stats -T \$fname_minlen | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv
        fi
        fname_subsampled=seqs_\${n_reads}_minlen_\${min_readlen}.fastq
        if [ ! -f \$fname_subsampled ]; then
            reformat.sh qin=33 ow=t samplereadstarget=\$n_reads in=\$fname_minlen out=\$fname_subsampled
            echo -n "subsample_\${n_reads}_minlen_\${min_readlen}\t" >> stats.tsv
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
            corOutCoverage=$cor_out_coverage \
            maxThreads=$task.cpus \
            maxMemory=${task.memory.toGiga()}g

        canu_status=\$?
        set -e
        if [[ \$canu_status -eq 0 ]] && [[ -f \$fname_contigs ]] && [[ -s \$fname_contigs ]]; then
            mv \$fname_contigs assemblies.fasta
            break
        else
            continue
        fi
    done

    # remove potentially large fastq files.
    rm minlen_*.fastq

    # creates empty file if not exist
    touch assemblies.fasta
    mv stats.tsv assembly_stats_${meta.mapping_target}.tsv
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


process canu_contig_info {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("contigs.fasta")
    output:
        tuple val(meta), path("contig-info-${meta.mapping_target}.tsv")
    """
    #!/usr/bin/env python
    import pandas as pd
    with open('contigs.fasta') as f_in:
        fasta_headers = [
            line[1:].strip().split()
            for line in f_in
            if line.startswith('>')
        ]
    info = [
        {
            'WorkflowMappingTarget': '${meta.mapping_target}',
            'Contig': header[0],
            **dict(entry.split('=') for entry in header[1:])
        }
        for header in fasta_headers
    ]
    pd.DataFrame(info).to_csv('contig-info-${meta.mapping_target}.tsv', sep='\\t', index=False)
    """
}


// process contigs_stats {
//     label "opr_general"
//     cpus 1
//     input:
//         tuple val(meta), path("sorted-contigs.fasta")
//     output:
//         tuple val(meta), path("stats.tsv")
//     """
//     seqkit stats sorted-contigs.fasta > stats.tsv
//     """
// }


// process contigs_readlengths {
//     label "opr_general"
//     cpus 1
//     input:
//         tuple val(meta), path("sorted-contigs.fasta")
//     output:
//         tuple val(meta), path("readlengths.txt")
//     """
//     readlength.sh qin=33 in=sorted-contigs.fasta bin=1 out=readlengths.txt
//     """
// }


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
        -outfmt 5 \
        -max_target_seqs 1
    """
}


process extract_blasthits {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta), path("blast-results.xml")
    output:
        tuple val(meta), path("blast-hits-${meta.mapping_target}.csv")
    """
    #!/usr/bin/env python
    import pandas as pd
    from xml.etree import ElementTree
    rows = []
    try:
        root = ElementTree.parse('blast-results.xml').getroot()
        for iteration in root.findall(".//Iteration"):
            query_name = iteration.find('Iteration_query-def').text.split()[0]
            hit = iteration.find("./Iteration_hits/Hit[Hit_num='1']")
            if hit is not None:
                try:
                    accession = hit.find('Hit_accession').text
                    fasta_header = hit.find('Hit_def').text.replace(',', '').replace(';', '')
                    # sometimes, the description contains the | symbol, so we have to split around it.
                    id_and_description, family, organism, orientation, segment = fasta_header.rsplit('|', 4)
                    description = id_and_description.split('|', 1)[1]
                    hit_length = int(hit.find('Hit_len').text)
                    hsp = hit.find("./Hit_hsps/Hsp[Hsp_num='1']")
                    if hsp is not None:
                        bit_score = float(hsp.find("Hsp_bit-score").text)
                        query_from = int(hsp.find("Hsp_query-from").text)
                        query_to = int(hsp.find("Hsp_query-to").text)
                        hit_from = int(hsp.find("Hsp_hit-from").text)
                        hit_to = int(hsp.find("Hsp_hit-to").text)
                        identity = int(hsp.find("Hsp_identity").text)
                        alignment_length = int(hsp.find("Hsp_align-len").text)
                        gaps = int(hsp.find("Hsp_gaps").text)
                        rows.append({
                            'Query': query_name.strip(),
                            'Reference': accession.strip(),
                            'Description': description.strip(),
                            'Family': family.strip(),
                            'Organism': organism.strip(),
                            'Segment': segment.strip(),
                            'Orientation': orientation.strip(),
                            'HitLength': hit_length,
                            'Bitscore': bit_score,
                            'QueryFrom': query_from,
                            'QueryTo': query_to,
                            'HitFrom': hit_from,
                            'HitTo': hit_to,
                            'IdenticalPositions': identity,
                            'AlignmentLength': alignment_length,
                            'Gaps': gaps,
                            'WorkflowMappingTarget': '${meta.mapping_target}',
                        })
                except (ValueError, AttributeError):
                    continue
    except ElementTree.ParseError:
        # empty file with no blast hits
        pass
    hits = pd.DataFrame(rows)
    hits.to_csv('blast-hits-${meta.mapping_target}.csv', header=True, index=False)
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
    minimap2 -ax map-ont ref.fasta trimmed.fastq -t ${task.cpus} | samtools view -F 4 -b -o mapped.bam
    samtools sort --threads ${task.cpus} -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}


process calc_coverage {
    label "opr_general"
    cpus 8
    input:
        tuple val(meta), path("sorted.bam"), path("sorted.bam.bai")
    output:
        tuple val(meta), path("coverage.txt")
    """
    samtools depth -a sorted.bam > coverage.txt
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
            path(model),
            val(min_depth)
    output:
        tuple val(meta), path('consensus.fasta')
    """
    medaka inference sorted.bam consensus.hdf --threads ${task.cpus} --model ${model}
    medaka sequence consensus.hdf ref.fasta consensus.fasta --fill_char N --min_depth ${min_depth}
    """
}


process simple_consensus {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            val(min_depth),
            val(min_share)
    output:
        tuple val(meta), path("cons.fasta")
    """
    samtools consensus \
    -a --mark-ins --show-del yes --show-ins yes \
    --mode simple \
    --min-depth $min_depth \
    --call-fract $min_share \
    -f fasta sorted.bam \
    > consensus_draft.fasta

    # write a consensus of Ns in case there was nothing mapped at all.
    if [ ! -s consensus_draft.fasta ]; then
        header=\$(grep ">" ref.fasta)
        seq_length=\$(readlength.sh ref.fasta | grep -w '#Bases:' | awk '{print \$2}')
        sequence=\$(printf "%.0sN" \$(seq 1 "\$seq_length") | fold -w 60)
        echo \$header > consensus_draft.fasta
        echo \$sequence >> consensus_draft.fasta
    fi

    sed 's/^> *\\([^ ]*\\)/>Mapped_to_\\1/' consensus_draft.fasta > consensus_draft_header_corrected.fasta

    consensus_correction ref.fasta consensus_draft_header_corrected.fasta sorted.bam \
    --call-fract $min_share \
    --min-depth $min_depth \
    --output cons.fasta
    """
}


process compute_mapping_stats {
    label "opr_general"
    cpus 1
    input:
        tuple val(meta),
        path('ref.fasta'),
        path('sorted.bam'),
        path('sorted.bam.bai'),
        path('cons.fasta'),
        path('depth.txt')
    output:
        tuple val(meta), path("mapping_stats.csv")
    """
    ref_id=\$(head -n 1 ref.fasta | sed 's/>//' | awk '{print \$1}')
    ref_len=\$(seqtk size ref.fasta | awk '{print \$2}')
    num_mapped_reads=\$(samtools view -F 4 -c sorted.bam)
    n_count=\$(seqtk comp cons.fasta | awk '{print \$9}')
    cons_length=\$(seqtk comp cons.fasta | awk '{print \$2}')
    avg_coverage=\$(awk '{sum+=\$3} END {if (NR > 0) print sum/NR; else print 0}' depth.txt)
    echo "\$ref_id\t\$ref_len\t\$num_mapped_reads\t\$n_count\t\$cons_length\t\$avg_coverage" > mapping_stats.csv
    """
}


process concat_mapping_stats {
    label "opr_general"
    cpus 1
    input:
        tuple val(samplename), path("collected_stats_*.tsv")
    output:
        tuple val(samplename), path('all_stats.tsv')
    """
    echo "Reference\tReferenceLength\tNumberOfMappedReads\tNCount\tConsensusLength\tAverageCoverage" > all_stats.tsv
    cat collected_stats_*.tsv >> all_stats.tsv
    """
}


process sample_report {
    label "wf_common"
    cpus 1
    input:
        tuple val(samplename),
            path('clean_stats.tsv'),
            val(assembly_modes),
            path(assembly_stats),
            path(blast_hits),
            path('mapping_stats.tsv'),
            path(contig_infos),
            path('virus_db_config.yaml')
    output:
        tuple val(samplename), path("${samplename}.html"), path("${samplename}_consensus_stats.tsv")
    """
    workflow-glue report_sample \
        --prefix ${samplename} \
        --outdir . \
        --mapping-stats mapping_stats.tsv \
        --clean-read-stats clean_stats.tsv \
        --assembly-modes ${assembly_modes.join(" ")} \
        --assembly-read-stats ${assembly_stats.join(" ")} \
        --blast-hits ${blast_hits.join(" ")} \
        --contig-info ${contig_infos.join(" ")} \
        --virus-db-config virus_db_config.yaml
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
        trimmed = samples | map{ meta, reads, stats -> [meta, reads] } | trim
        trim_stats = trimmed | seqtk_stats
        trim_fastqc = trimmed | fastqc

        // metagenomic read classification with centrifuge
        classification = trimmed
        | combine(Channel.of(params.classification_db))
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
        assemblies = to_assemble_notarget
        | mix(to_assemble_targeted)
        | assemble_canu

        // database search for references using blast
        blast_queries = assemblies
        | map { meta, contigs, stats -> [meta, contigs] }
        | prepare_blast_search

        collected_contigs_infos = blast_queries
        | canu_contig_info
        | map { meta, contig_info -> [meta.alias, contig_info]}
        | groupTuple(by: 0)

        blast_hits = blast_queries
        | map { meta, contigs -> [meta, contigs, params.blast_db, params.all_target] }
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
        | map { ref_id -> [ref_id, params.blast_db, params.all_target]}
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
        mapped_to_ref = reads_and_ref | map_to_ref

        // get coverages
        coverage = mapped_to_ref | map { meta, ref, bam, bai -> [meta, bam, bai]} | calc_coverage

        if (params.consensus_method == 'simple') {
            consensi = mapped_to_ref
            | map { meta, ref, bam, bai -> [meta, ref, bam, bai, params.consensus_min_depth, params.consensus_min_share] }
            | simple_consensus
        } else if (params.consensus_method == 'medaka') {
            consensi = mapped_to_ref
            | map { meta, ref, bam, bai -> [meta, ref, bam, bai, params.medaka_consensus_model, params.consensus_min_depth] }
            | medaka_consensus
        }

        // TODO: create & export vcf files?
        // TODO: replace header in consensus fasta files! (Do already in consensus step!)

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

        // Create the report.
        collected_clean_stats = cleaned
        | map { meta, reads, stats -> [meta.alias, stats] }

        assembly_modes = params.assemble_notarget ? params.targets + ['no-target'] : params.targets

        collected_assembly_stats = assemblies
        | map {meta, contigs, stats -> [meta.alias, stats]}
        | groupTuple(by: 0)

        collected_blast_hits = blast_hits
        | map {meta, hits -> [meta.alias, hits]}
        | groupTuple(by: 0)

        sample_reports = collected_clean_stats
        | map {samplename, clean_stats -> [samplename, clean_stats, assembly_modes]}
        | join(collected_assembly_stats, by: 0)
        | join(collected_blast_hits, by: 0)
        | join(collected_mapping_stats, by: 0)
        | join(collected_contigs_infos, by: 0)
        | combine(Channel.of(params.virus_db_config))
        | sample_report

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
            mapped_to_ref | map { meta, ref, bam, bai -> [ref, "$meta.alias/consensus/mapping", "${meta.consensus_target}.fasta"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bam, "$meta.alias/consensus/mapping", "${meta.consensus_target}.bam"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bai, "$meta.alias/consensus/mapping", "${meta.consensus_target}.bam.bai"] },
            consensi | map { meta, consensus -> [consensus, "$meta.alias/consensus/consensus/${meta.consensus_target}", null] },
            coverage | map { meta, coverage -> [coverage, "$meta.alias/consensus/consensus/${meta.consensus_target}", null] },
            // Report related output
            collected_mapping_stats | map { alias, stats -> [stats, "$alias/sample_stats", null]},
            sample_reports | map { alias, html_report, consensus_stats -> [html_report, "$alias/report", null]},
            sample_reports | map { alias, html_report, consensus_stats -> [consensus_stats, "$alias/report", null]},
            // advanced output
            cleaned | filter { params.advanced_output.cleaned_reads } | map { meta, reads, stats -> [reads, "$meta.alias/clean", null] },
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
