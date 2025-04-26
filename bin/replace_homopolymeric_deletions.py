#!/usr/bin/env python3

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
            counter = 0
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


def process_vcf_set_to_N(
        in_vcf_path,
        out_vcf_path,
        max_del_len,
        homopolymeric_indices,
        remove
):
    '''
    Read variants from in_vcf_path, replace qualifying deletions in homopolymers
    with N-substitutions, and write to out_vcf_path.
    '''
    with pysam.VariantFile(in_vcf_path, 'r') as vcf_in:
        with pysam.VariantFile(out_vcf_path, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                new_alts = []
                if record.alts:
                    for alt in record.alts:
                        if is_homopolymeric_deletion(record, alt, max_del_len, homopolymeric_indices):
                            if not remove:
                                # Build a substitution: keep the first (anchor) base, then Ns
                                anchor = record.ref[0]
                                del_len = len(record.ref) - len(alt)
                                new_alt = anchor + ('N' * del_len)
                                new_alts.append(new_alt)
                        else:
                            new_alts.append(alt)
                if new_alts:
                    # Replace the ALT alleles in the record
                    record.alts = tuple(new_alts)
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
    parser.add_argument(
        '--remove',
        required=True,
        help='Whether to remove the deletions from the homopolymeric set or replace them with and N (yes for removal)'
    )
    args = parser.parse_args()

    ref = SeqIO.read(args.reference, 'fasta')
    homopolymeric_regions = get_homopolymeric_regions(
        ref.seq, args.min_homopolymer_length
    )

    process_vcf_set_to_N(
        in_vcf_path=args.vcf,
        out_vcf_path=args.out,
        max_del_len=args.max_deletion_length,
        homopolymeric_indices=homopolymeric_regions,
        remove=(args.remove=='yes')
    )


if __name__ == '__main__':
    main()
