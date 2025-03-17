nextflow.enable.dsl = 2


process lengths_and_qualities {
    label "general"
    cpus 1
    input:
        tuple val(samplename), path(reads)
    output:
        tuple val(samplename), path("read_len_qual.tsv")
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from pandas import DataFrame
    lines = [
        (len(record.seq), sum(record.letter_annotations['phred_quality']) / len(record.seq))
        for record in SeqIO.parse('${reads}', 'fastq')
        if len(record.seq) > 0
    ]
    DataFrame(lines, columns=['Length', 'Quality']).to_csv('read_len_qual.tsv', sep='\t')
    """
}


process trim {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("demultiplexed.fastq.gz")
    output:
        tuple val(meta), path("trimmed.fastq")
    """
    seqtk trimfq -b ${params.trim_length} -e ${params.trim_length} demultiplexed.fastq.gz > trimmed.fastq
    """
}


process classify_centrifuge {
    label "centrifuge"
    cpus 12
    memory '28 GB'
    input:
        tuple val(meta), path("seqs.fastq"), path(db_path), val(target_db)
    output:
        tuple val(meta),
            path("classification_${target_db}.tsv"),
            path("classification_report_${target_db}.tsv"),
            path("classification_kraken_${target_db}.tsv"),
            path("classification_${target_db}.html")
    """

    if [ -s seqs.fastq ]
    then
        centrifuge \
            -p ${task.cpus} \
            --mm \
            -x ${db_path}/${target_db} \
            -U seqs.fastq \
            --report-file classification_report_${target_db}.tsv \
            -S classification_${target_db}.tsv
        centrifuge-kreport \
            -x ${db_path}/${target_db} classification_${target_db}.tsv > classification_kraken_${target_db}.tsv
        ktImportTaxonomy \
            -tax ${db_path}/taxonomy \
            -m 3 -t 5 \
            classification_kraken_${target_db}.tsv \
            -o classification_${target_db}.html
    else
        # write a report even if no sequences were available
        touch classification_${target_db}.tsv
        touch classification_report_${target_db}.tsv
        touch classification_kraken_${target_db}.tsv

        # Create an informative HTML file telling the user that no sequences were found
        cat <<EOF > classification_${target_db}.html
<!DOCTYPE html>
<html>
<head>
    <title>No Sequences Found</title>
    <style>
        body { font-family: Arial, sans-serif; text-align: center; margin-top: 50px; }
        h1 { color: red; }
    </style>
</head>
<body>
    <h1>No sequences were found</h1>
    <p>For this barcode no sequences could be classified by centrifuge.</p>
</body>
</html>
EOF
    fi
    """
}


process filter_contaminants {
    label "general"
    cpus 4
    memory '20 GB'
    input:
        tuple val(meta), path('seqs.fastq'), path('db_*.fna.gz'), val(contaminants)
    output:
        tuple val(meta), path('filtered.fastq'), emit: reads
        tuple val(meta.alias), path("clean_stats.tsv"), emit: stats
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
    label "general"
    cpus 4
    memory '10 GB'
    input:
        tuple val(meta), path('seqs.fastq'), path(target)
    output:
        tuple val(meta), path('filtered.fastq')
    """
    minimap2 \
        -x map-ont \
        -a ${target} \
        -t ${task.cpus} \
        seqs.fastq \
    | samtools fastq \
        -F 4 > filtered.fastq
    """
}

process assemble_canu {
    label "canu"
    cpus 16
    memory '24 GB'
    input:
        tuple val(meta),
            path("seqs.fastq"),
            val(min_read_length),
            val(min_overlap_length),
            val(read_sampling_bias),
            val(genome_size),
            val(cor_out_coverage),
            val(stop_on_low_coverage),
            val(min_input_coverage),
            val(max_input_coverage)
    output:
        tuple val(meta), path("asm.contigs.fasta"), emit: contigs
        tuple val(meta), path("assembly_stats_${meta.mapping_target}.tsv"), emit: stats
    """        
    echo "step\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > stats.tsv
    echo -n "all_${meta.mapping_target}\t" >> stats.tsv
    if [ -s seqs.fastq ]
    then
        seqkit stats -T seqs.fastq | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv
    else
        echo "0\t0\t0\t0.0\t0" >> stats.tsv
    fi

    outdir=.

    set +e

    canu \
        -nanopore-raw seqs.fastq \
        -fast \
        -p asm \
        -d \$outdir \
        genomeSize=$genome_size \
        minReadLength=$min_read_length  \
        minOverlapLength=$min_overlap_length \
        corOutCoverage=$cor_out_coverage \
        readSamplingBias=$read_sampling_bias \
        stopOnLowCoverage=$stop_on_low_coverage \
        minInputCoverage=$min_input_coverage \
        maxInputCoverage=$max_input_coverage \
        maxThreads=${task.cpus} \
        maxMemory=${task.memory.toGiga()}g

    canu_status=\$?  # Capture Canu's exit code

    set -e

    touch asm.contigs.fasta

    if [[ ! -f asm.correctedReads.fasta.gz ]]
    then
        touch asm.correctedReads.fasta
        gzip asm.correctedReads.fasta
    fi

    echo -n "corrected_${meta.mapping_target}\t" >> stats.tsv

    first_line_count=\$(gunzip -c asm.correctedReads.fasta.gz | head -n 1 | wc | awk '{print \$1}')

    if [[ \$first_line_count -ne 0 ]]
    then
        seqkit stats -T asm.correctedReads.fasta.gz | tail -n 1 | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> stats.tsv
    else
        echo "0\t0\t0\t0.0\t0" >> stats.tsv
    fi

    mv stats.tsv assembly_stats_${meta.mapping_target}.tsv
    """
}


process pop_bubbles {
    label "general"
    cpus 1
    input:
        tuple val(meta), path('canu_contigs.fasta')
    output:
        tuple val(meta), path('nobubbles.fasta')
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    with open('nobubbles.fasta', "w") as f_out:
        for record in SeqIO.parse('canu_contigs.fasta', 'fasta'):
            if "bubble=yes" not in record.description:
                SeqIO.write(record, f_out, 'fasta')
    """
}


process prepare_blast_search {
    label "general"
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
    label "general"
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
    columns = [
        'WorkflowMappingTarget',
        'Contig',
        'len',
        'reads',
    ]
    df = pd.DataFrame(info, columns=columns)
    df.to_csv('contig-info-${meta.mapping_target}.tsv', sep='\\t', index=False)
    """
}


process blast {
    label "general"
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
    label "general"
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
        hits = pd.DataFrame(rows)
    except ElementTree.ParseError:
        # empty file with no blast hits
        column_names = [
            'Query',
            'Reference',
            'Description',
            'Family',
            'Organism',
            'Segment',
            'Orientation',
            'HitLength',
            'Bitscore',
            'QueryFrom',
            'QueryTo',
            'HitFrom',
            'HitTo',
            'IdenticalPositions',
            'AlignmentLength',
            'Gaps',
            'WorkflowMappingTarget',
        ]
        hits = pd.DataFrame(columns=column_names)
    hits.to_csv('blast-hits-${meta.mapping_target}.csv', header=True, index=False)
    """
}


process get_ref_fasta {
    label "general"
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
    label "general"
    cpus 8
    memory '10 GB'
    input:
        tuple val(meta), path("trimmed.fastq"), path("ref.fasta")
    output:
        tuple val(meta), path("ref.fasta"), path("sorted.bam"), path("sorted.bam.bai")
    """
    minimap2 -ax map-ont ref.fasta trimmed.fastq \
        -t ${task.cpus} --secondary=no \
        -w ${params.map_to_target_minimap_window_size} \
        -k ${params.map_to_target_minimap_kmer_size} \
    | samtools view -F 4 -b -o mapped.bam

    samtools sort --threads ${task.cpus} -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}


process calc_coverage {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("sorted.bam"), path("sorted.bam.bai")
    output:
        tuple val(meta), path("coverage.txt")
    """
    samtools depth -aa -J sorted.bam > coverage.txt
    """
}


process simplify_reference_fasta {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("ref.fasta")
    output:
        tuple val(meta), path("simple_ref.fasta")
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    ref = SeqIO.read('ref.fasta', 'fasta')
    with open('simple_ref.fasta', 'w') as f_out:
        f_out.write(f">{ref.id}\\n{ref.seq}\\n") 
    """
}


process medaka_variant_consensus {
    label "medaka"
    cpus 2
    memory '24 GB'
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            path("trimmed.fastq"),
            val(model),
            val(min_depth)
    output:
        tuple val(meta), path("consensus.${meta.consensus_target}.fasta"), emit: consensus
        tuple val(meta), path("variants.${meta.consensus_target}.vcf.gz"), emit: variants
        tuple val(meta), path("depth.${meta.consensus_target}.txt"), emit: depth
    """
    samtools depth -aa -J sorted.bam > depth.txt

    medaka_cmd="medaka inference sorted.bam consensus.hdf --threads ${task.cpus}"

    model_choice="default"

    if [[ ${model} == "auto" ]]
    then
        set +e
        model_path=\$(medaka tools resolve_model --auto_model variant trimmed.fastq 2>/dev/null)
        exit_code=\$?
        if [[ \$exit_code -eq 0 && -n "\$model_path" ]]
        then
            medaka_cmd="\$medaka_cmd --model \$model_path"
            model_choice=\$model_path
        fi
        set -e
    else
        medaka_cmd="\$medaka_cmd --model ${model}:variant"
        model_choice=${model}:variant
    fi

    \$medaka_cmd

    medaka vcf consensus.hdf ref.fasta variants.vcf

    bcftools sort variants.vcf -o sorted.vcf
    medaka tools annotate sorted.vcf ref.fasta sorted.bam annotated.vcf

    bcftools filter -e "ALT='.'" annotated.vcf \\
    | bcftools filter -o filtered.vcf -O v -e "INFO/DP<${min_depth}" -

    refid=\$(head -n 1 ref.fasta | awk '{print substr(\$1, 2)}')

    awk -v ref="\$refid" -v thresh=${min_depth} '{if (\$3<thresh) print \$1"\\t"\$2+1}' depth.txt  > mask.regions

    bgzip filtered.vcf
    tabix filtered.vcf.gz

    bcftools consensus \\
        --mask mask.regions \\
        --fasta-ref ref.fasta \\
        -o consensus_draft.fasta filtered.vcf.gz

    echo ">consensus method=medaka reference=\$refid sample=${meta.alias} medaka_model=\$model_choice" > consensus.fasta

    # write a consensus of Ns in case there was nothing mapped at all.
    if [ ! -s consensus_draft.fasta ]
    then
        seq_length=\$(readlength.sh ref.fasta | grep -w '#Bases:' | awk '{print \$2}')
        sequence=\$(printf "%.0sN" \$(seq 1 "\$seq_length") | fold -w 80)
        echo "\$sequence" >> consensus.fasta
    else
        seqkit seq -w 80 consensus_draft.fasta | tail -n +2 >> consensus.fasta
    fi

    mv depth.txt depth.${meta.consensus_target}.txt
    mv consensus.fasta consensus.${meta.consensus_target}.fasta
    mv filtered.vcf.gz variants.${meta.consensus_target}.vcf.gz
    """
}


process medaka_consensus {
    label "medaka"
    cpus 2
    memory '24 GB'
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            path("trimmed.fastq"),
            val(model),
            val(min_depth)
    output:
        tuple val(meta), path("consensus.${meta.consensus_target}.fasta")
    """
    medaka_cmd="medaka inference sorted.bam consensus.hdf --threads ${task.cpus}"

    model_choice="default"

    if [[ ${model} == "auto" ]]
    then
        set +e
        model_path=\$(medaka tools resolve_model --auto_model consensus trimmed.fastq 2>/dev/null)
        exit_code=\$?
        if [[ \$exit_code -eq 0 && -n "\$model_path" ]]
        then
            medaka_cmd="\$medaka_cmd --model \$model_path"
            model_choice=\$model_path
        fi
        set -e
    else
        medaka_cmd="\$medaka_cmd --model ${model}:consensus"
        model_choice=${model}:consensus
    fi

    \$medaka_cmd

    medaka sequence consensus.hdf ref.fasta consensus_draft.fasta --fill_char N --min_depth ${min_depth}

    refid=\$(head -n 1 ref.fasta | awk '{print substr(\$1, 2)}')
    echo ">consensus method=medaka reference=\$refid sample=${meta.alias} medaka_model=\$model_choice" > consensus.fasta

    # write a consensus of Ns in case there was nothing mapped at all.
    if [ ! -s consensus_draft.fasta ]
    then
        seq_length=\$(readlength.sh ref.fasta | grep -w '#Bases:' | awk '{print \$2}')
        sequence=\$(printf "%.0sN" \$(seq 1 "\$seq_length") | fold -w 80)
        echo "\$sequence" >> consensus.fasta
    else
        seqkit seq -w 80 consensus_draft.fasta | tail -n +2 >> consensus.fasta
    fi

    mv consensus.fasta consensus.${meta.consensus_target}.fasta
    """
}


process simple_consensus {
    label "general"
    cpus 1
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            val(min_depth),
            val(min_share)
    output:
        tuple val(meta), path("consensus.${meta.consensus_target}.fasta")
    """
    samtools consensus \
    -a --mark-ins --show-del yes --show-ins yes \
    --mode simple \
    --min-depth $min_depth \
    --call-fract $min_share \
    -f fasta sorted.bam \
    > consensus_draft.fasta

    refid=\$(head -n 1 ref.fasta | awk '{print substr(\$1, 2)}')
    echo ">consensus method=samtools_simple reference=\$refid sample=${meta.alias}" > consensus.fasta

    # write a consensus of Ns in case there was nothing mapped at all.
    if [ ! -s consensus_draft.fasta ]
    then
        seq_length=\$(readlength.sh ref.fasta | grep -w '#Bases:' | awk '{print \$2}')
        sequence=\$(printf "%.0sN" \$(seq 1 "\$seq_length") | fold -w 80)
        echo "\$sequence" >> consensus.fasta
    else
        seqkit seq -w 80 consensus_draft.fasta | tail -n +2 >> consensus.fasta
    fi

    correct_samtools_consensus.py ref.fasta consensus.fasta sorted.bam \
    --call-fract $min_share \
    --min-depth $min_depth \
    --output consensus.${meta.consensus_target}.fasta
    """
}


process compute_mapping_stats {
    label "general"
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


process empty_tsv {
    label "general"
    cpus 1
    input:
        val(samplename)
    output:
        tuple val(samplename), path('empty.tsv')
    """
    touch empty.tsv
    """
}


process concat_mapping_stats {
    label "general"
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


process empty_fasta {
    label "general"
    cpus 1
    output:
        path('empty.fasta')
    """
    touch empty.fasta
    """
}


process collect_reference_info {
    label "general"
    cpus 1
    input:
        path('ref_*.fasta')
    output:
        path('reference_info.tsv')
    """
    #!/usr/bin/env python
    import pandas as pd
    import glob
    data = []
    empty = ''
    for fname in glob.glob("ref_*.fasta"):
        with open(fname) as f:
            header = f.readline().lstrip('>').strip()
            if not header:
                continue
            try:
                id, rest = map(str.strip, header.split('|', 1))
                descr, fam, org, orient, seg = map(str.strip, rest.rsplit('|', 4))
                data.append({
                    'Reference': id,
                    'Description': descr,
                    'Family': fam,
                    'Organism': org,
                    'Orientation': orient,
                    'Segment': seg,
                })
            except ValueError:
                data.append({
                    'Reference': header.split('|')[0],
                    'Description': empty,
                    'Family': empty,
                    'Organism': empty,
                    'Orientation': empty,
                    'Segment': empty,
                })
    cols = [
        'Reference',
        'Description',
        'Family',
        'Organism',
        'Orientation',
        'Segment',
    ]
    df = pd.DataFrame(data, columns=cols)
    df.to_csv('reference_info.tsv', sep='\t')
    """
}


process sample_report {
    label "report"
    cpus 1
    input:
        tuple val(samplename),
            path('clean_stats.tsv'),
            val(assembly_modes),
            path(assembly_stats),
            path(blast_hits),
            path('mapping_stats.tsv'),
            path(contig_infos),
            path('trimmed_read_stats.tsv'),
            path('cleaned_read_stats.tsv'),
            path('virus_db_config.yaml'),
            path('reference_info.tsv')
    output:
        tuple val(samplename), path('report.html'), emit: report
        tuple val(samplename), path('stats_reads.tsv'), emit: read_stats
        tuple val(samplename), path('stats_contigs.tsv'), emit: contig_stats
        tuple val(samplename), path('stats_consensus.tsv'), emit: consensus_stats
    """
    workflow-glue mergestats_reads \
        --clean-read-stats clean_stats.tsv \
        --out stats_reads.tsv

    workflow-glue mergestats_contig \
        --blast-hits ${blast_hits.join(" ")} \
        --contig-info ${contig_infos.join(" ")} \
        --out stats_contigs.tsv

    workflow-glue mergestats_consensus \
        --virus-db-config virus_db_config.yaml \
        --mapping-stats mapping_stats.tsv \
        --reference-info reference_info.tsv \
        --out stats_consensus.tsv

    workflow-glue sample_html_report \
        --pipeline-version ${workflow.manifest.version} \
        --samplename ${samplename} \
        --virus-db-config virus_db_config.yaml \
        --reads-stats stats_reads.tsv \
        --contigs-stats stats_contigs.tsv \
        --consensus-stats stats_consensus.tsv \
        --trimmed-read-distribution trimmed_read_stats.tsv \
        --cleaned-read-distribution cleaned_read_stats.tsv \
        --out report.html
    """
}


process get_best_consensus_files {
    label "report"
    cpus 1
    input:
        tuple val(samplename), path('consensus_stats.tsv'), path('cons_*.fasta')
    output:
        tuple val(samplename), path('out')
    """
    workflow-glue get_curated_consensus_genomes \
        --consensus-stats consensus_stats.tsv \
        --consensus-files cons_*.fasta \
        --out-dir out
    """
}


process output {
    label "general"
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
