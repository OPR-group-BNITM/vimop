#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.


import argparse
from Bio import SeqIO
import pysam


def get_homopolymeric_regions(seq, minlen):
    homopolymeric_positions = []
    counter = 0
    for i, nt in enumerate(seq[1:], 1):
        if nt == seq[i-1]:
            counter += 1
        else:
            counter = 1
        if counter == minlen:
            homopolymeric_positions.extend(range(i - minlen + 1, i + 1))
        elif counter > minlen:
            homopolymeric_positions.append(i)
    return set(homopolymeric_positions)



def is_homopolymeric_deletion(record, alt, max_del_len, homopolymeric_indices):
    '''
    Check if the given ALT allele represents a deletion in a homopolymeric region
    that meets the size thresholds.
    '''
    ref = record.ref
    start = record.start

    if len(ref) <= 1 or len(alt) != 1:
        # it's not a deletion
        return False

    if alt[0] == ref[0]:
        ref = ref[1:]
        start += 1
    elif alt[-1] == ref[-1]:
        ref = ref[:-1]
    else:
        # not a deletion
        return False

    if len(set(ref)) > 1:
        # the deleted part is not homopolymeric
        return False

    # Check length bounds
    if len(ref) > max_del_len:
        return False

    return start in homopolymeric_indices


def remmove_homopolymeric_dels(
        in_vcf_path,
        out_vcf_path,
        max_del_len,
        homopolymeric_indices,
):
    '''
    Read variants from in_vcf_path, replace qualifying deletions in homopolymers
    with N-substitutions, and write to out_vcf_path.
    '''
    with pysam.VariantFile(in_vcf_path, 'r') as vcf_in:
        with pysam.VariantFile(out_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                is_homopolymeric_del = (
                    bool(record.alts) and is_homopolymeric_deletion(record, record.alts[0], max_del_len, homopolymeric_indices)
                )
                if not is_homopolymeric_del:
                    vcf_out.write(record)


def main():
    parser = argparse.ArgumentParser(
        description='Replace deletions in homopolymeric regions with N substitutions'
    )
    parser.add_argument(
        '--vcf',
        required=True,
        help='Path to input VCF file'
    )
    parser.add_argument(
        '--reference',
        required=True,
        help='Path to reference genome (fasta)'
    )
    parser.add_argument(
        '-o', '--out',
        required=True,
        help='Path to output VCF file'
    )
    parser.add_argument(
        '--min-homopolymer-length',
        type=int,
        default=1,
        help='Minimum length of deletion to consider (homopolymer length threshold)'
    )
    parser.add_argument(
        '--max-deletion-length',
        type=int,
        default=10,
        help='Maximum deletion length to process'
    )
    args = parser.parse_args()

    ref = SeqIO.read(args.reference, 'fasta')
    homopolymeric_regions = get_homopolymeric_regions(
        ref.seq, args.min_homopolymer_length
    )

    remmove_homopolymeric_dels(
        in_vcf_path=args.vcf,
        out_vcf_path=args.out,
        max_del_len=args.max_deletion_length,
        homopolymeric_indices=homopolymeric_regions,
    )


if __name__ == '__main__':
    main()
