#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
import pandas as pd
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        description='Filter non-viral reads from Centrifuge output'
    )
    parser.add_argument(
        '--centrifuge',
        required=True,
        help='Centrifuge output file'
    )
    parser.add_argument(
        '--taxids',
        required=True,
        help='File with virus taxIDs (one per line)'
    )
    parser.add_argument(
        '--fastq',
        required=True,
        help='Original FASTQ file'
    )
    parser.add_argument(
        '--min-score',
        type=float,
        default=100.0,
        help='Minimum score to consider non-virus classification high-confidence (default: 100)'
    )
    parser.add_argument(
        '--out',
        default='filtered_reads.fastq',
        help='Output FASTQ file'
    )
    args = parser.parse_args()

    with open(args.taxids) as f:
        # 0 is unmapped
        taxids_virus = set(map(int, f))

    # get reads to remove
    df = pd.read_csv(args.centrifuge, sep='\t')
    mask = ((df['score'] > args.min_score) & (~df['taxID'].isin(taxids_virus.union(0))))
    remove = set(df.loc[mask, 'readID'])

    # filter the reads
    with open(args.out, 'w') as out_f:
        for record in SeqIO.parse(args.fastq, 'fastq'):
            if record.id not in remove:
                SeqIO.write(record, out_f, 'fastq')


if __name__ == '__main__':
    main()
