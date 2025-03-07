import subprocess
import re
import pytest
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam
import os
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from correct_samtools_consensus import get_valid_deletions, correct_consensus_indels


def compute_consensus(fname_bam, fname_consensus, min_depth, call_fraction):
    cmd = [
        'samtools', 'consensus',
        '-a',
        '--mark-ins',
        '--show-del', 'yes',
        '--show-ins', 'yes',
        '--mode', 'simple',
        '--min-depth', str(min_depth),
        '--call-fract', str(call_fraction),
        '-f', 'fasta',
        '-o', str(fname_consensus),
        str(fname_bam)
    ]
    try:
        subprocess.run(cmd, shell=False, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        cmd_str = ' '.join(cmd)
        msg = (
            f"Command '{cmd_str}' failed with return code {e.returncode}\n"
            f"Standard Output:\n{e.stdout}\n", 
            f"Standard Error:\n{e.stderr}\n", 
        )
        raise RuntimeError(msg)


def cigarstring_to_tuples(cigar_string):
    # Mapping of CIGAR characters to their corresponding integer codes
    cigar_ops = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    tuples = [(cigar_ops[op], int(length)) for length, op in pattern.findall(cigar_string)]
    return tuples


def make_bam_file(bam_entries, fname, reflen):
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': reflen, 'SN': 'refseq'}]
    }
    output_bam = pysam.AlignmentFile(fname, "wb", header=header)
    bam_entries_sorted = sorted(bam_entries, key=lambda x: x['ref_start'])
    for i, bam in enumerate(bam_entries_sorted, 1):
        align = pysam.AlignedSegment()
        align.query_name = f'read_{i}'
        align.query_sequence = bam['seq']
        align.flag = bam['flag']
        align.cigar = cigarstring_to_tuples(bam['cigar'])
        align.reference_id = 0
        align.reference_start = bam['ref_start']
        output_bam.write(align)
    output_bam.close()


seq_a = 'GG'
seq_b = 'GCTGGGCTCA'
seq_c = 'TTTT'
seq_d = 'CATCAA'
seq_e = 'AG'

big_ref = seq_a + seq_b + seq_c + seq_d + seq_e
bam_big_ref_prefix = {
    'seq': seq_b,
    'flag': 0,
    'cigar': f'{len(seq_b)}=',
    'ref_start': len(seq_a),
}
bam_big_ref_suffix = {
    'seq': seq_d,
    'flag': 0,
    'cigar': f'{len(seq_d)}=',
    'ref_start': len(seq_a) + len(seq_b) + len(seq_c),
}
bam_big_ref_del = {
    'seq': seq_b + seq_d,
    'flag': 0,
    'cigar': f'{len(seq_b)}=' + f'{len(seq_c)}D' + f'{len(seq_d)}=',
    'ref_start': len(seq_a),
}

small_ref = seq_a + seq_b + seq_d + seq_e
bam_small_ref_prefix = {
    'seq': seq_b,
    'flag': 0,
    'cigar': f'{len(seq_b)}=',
    'ref_start': len(seq_a),
}
bam_small_ref_suffix = {
    'seq': seq_d,
    'flag': 0,
    'cigar': f'{len(seq_d)}=',
    'ref_start': len(seq_a) + len(seq_b)
}
bam_small_ref_ins = {
    'seq': seq_b + seq_c + seq_d,
    'flag': 0,
    'cigar': f'{len(seq_b)}=' + f'{len(seq_c)}I' + f'{len(seq_d)}=',
    'ref_start': len(seq_a),
}

testdata_integrate_map = {
    'no_deletion_called': {
        'refseq': big_ref,
        'bam_entries': (
            [bam_big_ref_del] * 10
            + [bam_big_ref_prefix] * 20
            + [bam_big_ref_suffix] * 20
        ),
        'min_depth': 20,
        'call_fract': 0.7,
        'expected_consensus': (len(seq_a) * 'N') + seq_b + (len(seq_c) * 'N') + seq_d + (len(seq_e) * 'N')
    },
    'deletion_called': {
        'refseq': big_ref,
        'bam_entries': (
            [bam_big_ref_del] * 10
            + [bam_big_ref_suffix] * 20
            + [bam_big_ref_del] * 20
        ),
        'min_depth': 5,
        'call_fract': 0.7,
        'expected_consensus': (len(seq_a) * 'N') + seq_b + seq_d + (len(seq_e) * 'N')
    },
    'no_insertion_called': {
        'refseq': small_ref,
        'bam_entries': (
            [bam_small_ref_ins] * 10
            + [bam_small_ref_prefix] * 20
            + [bam_small_ref_suffix] * 20
        ),
        'min_depth': 20,
        'call_fract': 0.7,
        'expected_consensus': (len(seq_a) * 'N') + seq_b + seq_d + (len(seq_e) * 'N')
    },
    'insertion_called': {
        'refseq': small_ref,
        'bam_entries': (
            [bam_small_ref_ins] * 10
            + [bam_small_ref_prefix] * 20
            + [bam_small_ref_suffix] * 20
        ),
        'min_depth': 5,
        'call_fract': 0.7,
        'expected_consensus': (len(seq_a) * 'N') + seq_b + seq_c + seq_d + (len(seq_e) * 'N')
    },
}


def _test_consensus_corrections(
        tmpdir,
        refseq,
        bam_entries,
        min_depth,
        call_fract,
        expected_consensus
):
    fname_alignments = tmpdir.join('alignments.bam')
    fname_reference = tmpdir.join('reference.fasta')
    fname_consensus = tmpdir.join('consensus.fasta')

    with open(fname_reference, "w") as f_out:
        ref_record = SeqRecord(Seq(refseq), id="refseq")
        SeqIO.write(ref_record, f_out, "fasta")

    make_bam_file(bam_entries, fname_alignments, len(refseq))
    compute_consensus(fname_alignments, fname_consensus, min_depth, call_fract)

    consensus = next(SeqIO.parse(fname_consensus, "fasta"))
    valid_dels = get_valid_deletions(str(fname_alignments), str(fname_reference), min_depth, call_fract)
    consensus_corrected_str, _, _ = correct_consensus_indels(str(consensus.seq), valid_dels, False)

    assert consensus_corrected_str == expected_consensus


@pytest.mark.parametrize("args", testdata_integrate_map.values(), ids=testdata_integrate_map.keys())
def test_consensus_corrections(tmpdir, args):
    _test_consensus_corrections(tmpdir, **args)
