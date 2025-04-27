#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

from collections import Counter
import argparse
import pysamstats
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def call(record, min_cov, min_frac):
    n_deletions = record['deletions']
    n_reads_all = record['reads_all']

    if n_reads_all < min_cov:
        return 'N'

    if n_deletions / n_reads_all >= min_frac:
        return ''

    acgt_counter = Counter({nt: record[nt] for nt in 'ACGT'})
    nt, count = acgt_counter.most_common(1)[0]
    if count / n_reads_all >= min_frac:
        return nt

    return 'N'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', required=True, help='Reference in fasta format')
    parser.add_argument('--bam', required=True, help='Bam file')
    parser.add_argument('--min-cov', required=True, type=int, help='Minimum coverage')
    parser.add_argument('--min-frac', required=True, type=float, help='Minimum fraction')
    parser.add_argument('--header', required=True, help='')
    parser.add_argument('--out', required=True, help='Fasta filename for output consensus')
    args = parser.parse_args()

    ref = SeqIO.read(args.ref, 'fasta')
    stats = {
        record['pos']: record
        for record in pysamstats.stat_variation(args.bam, fafile=args.ref)
    }
    consensus = ''.join(
        (call(stats[i], args.min_cov, args.min_frac) if i in stats else 'N')
        for i in range(len(ref.seq))
    )
    consensus_record = SeqRecord(
        id='consensus',
        name='consensus',
        description=args.header,
        seq=Seq(consensus)
    )
    SeqIO.write(
        consensus_record,
        args.out,
        'fasta'
    )    


if __name__ == '__main__':
    main()
