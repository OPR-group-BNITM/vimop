nextflow.enable.dsl = 2


process lengths_and_qualities {
    label "opr_general"
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
