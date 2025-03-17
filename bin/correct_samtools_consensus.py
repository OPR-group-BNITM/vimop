#!/usr/bin/env python3

import argparse
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysamstats


def get_valid_deletions(fname_bam, fname_reference, min_depth, min_fraction):
    """Get all deletions passing the minimum depth and fraction requirements.
    """
    return {
        record['pos']
        for record in pysamstats.stat_variation(fname_bam, fafile=fname_reference)
        if record['reads_all'] >= min_depth and record['deletions'] / record['reads_all'] >= min_fraction
    }


def correct_consensus_indels(consensus, valid_deletions, collect_indels):
    """Correct the consensus computed with samtools consensus.

    Samtools consensus also includes insertions and deletions below the minimum depth requirements.
    We do not want them included and instead keep a N.

    """
    # Remove all _N substrings from the sequence
    consensus = consensus.replace('_N', '')
    insertion_positions = [
        k - j
        for j, k in enumerate(i for i, c in enumerate(consensus) if c == '_')
    ]
    insertion_positions_set = set(insertion_positions)
    ref_positions = np.cumsum([
        1 if i not in insertion_positions_set else 0
        for i in range(len(consensus))
    ]) - 1
    # remove remaining insertion markers
    consensus = consensus.replace('_', '')
    deletion_positions = [
        i
        for (i, c), ref_pos in zip(enumerate(consensus), ref_positions)
        if c == '*' and int(ref_pos) in valid_deletions
    ]
    # replace the valid deletions with ''
    deletion_positions_set = set(deletion_positions)
    modified_consensus = ''.join([
        c
        for i, c in enumerate(consensus)
        if i not in deletion_positions_set
    ])
    # replace the invalid deletions with an N
    modified_consensus = modified_consensus.replace('*', 'N')

    if collect_indels and len(ref_positions) > 0:
        vcf_ins_candidates = [[] for _ in range(ref_positions[-1] + 1)]
        for ref_pos, nt in zip(ref_positions, consensus):
            vcf_ins_candidates[ref_pos].append(nt)
        vcf_ins = [
            (ref_pos, ''.join(nts))
            for ref_pos, nts in enumerate(vcf_ins_candidates)
            if len(nts) > 1
        ]
        vcf_del = [ref_positions[i] for i in deletion_positions]
    else:
        vcf_ins = []
        vcf_del = []

    return modified_consensus, vcf_ins, vcf_del


def create_vcf(reference_sequence, deletions, insertions, fname_vcf):
    with open(fname_vcf, 'w') as vcf_file:
        # Write the VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(f"##reference={reference_sequence.id}\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write the deletion variants
        for del_pos in deletions:
            ref_base = reference_sequence.seq[del_pos-1:del_pos+1]
            alt_base = reference_sequence.seq[del_pos-1]
            vcf_file.write(f"{reference_sequence.id}\t{del_pos}\t.\t{ref_base}\t{alt_base}\t.\t.\t.\n")

        # Write the insertion variants
        for ins_pos, ins_seq in insertions:
            ref_base = reference_sequence.seq[ins_pos]
            vcf_file.write(
                f"{reference_sequence.id}\t{ins_pos + 1}\t.\t{ref_base}\t{ins_seq}\t.\t.\t.\n"
            )


def main():
    parser = argparse.ArgumentParser(description="Process consensus and coverage based on a reference.")
    parser.add_argument("reference", help="FASTA file containing the reference sequence.")
    parser.add_argument("consensus", help="FASTA file containing the consensus sequence.")
    parser.add_argument("alignments", help="Alignment file (bam).")
    parser.add_argument("--min-depth", type=int, help="Minimum depth for considering deletions.", default=20)
    parser.add_argument("--call-fract", type=float, help="Minimum fraction for considering deletions.", default=0.6)
    parser.add_argument("--output", "-o", help="Output FASTA file for the modified consensus sequence.", default="consensus_corrected.fasta")
    parser.add_argument("--vcf")
    args = parser.parse_args()

    consensus = next(SeqIO.parse(args.consensus, "fasta"))
    valid_dels = get_valid_deletions(args.alignments, args.reference, args.min_depth, args.call_fract)
    consensus_corrected_str, vcf_ins, vcf_del = correct_consensus_indels(str(consensus.seq), valid_dels, bool(args.vcf))
    consensus_record = SeqRecord(
        Seq(consensus_corrected_str),
        id=consensus.id,
        description=consensus.description
    )
    with open(args.output, "w") as f_out:
        SeqIO.write(consensus_record, f_out, "fasta")
    if args.vcf:
        refseq = SeqIO.read(args.reference, "fasta")
        create_vcf(refseq, vcf_del, vcf_ins, args.vcf)


if __name__ == '__main__':
    main()
