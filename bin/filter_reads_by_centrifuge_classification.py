#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
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
        '--virus-taxids',
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

    with open(args.virus_taxids) as f:
        # 0 is unclassified
        taxids_virus = set(int(line.strip()) for line in f if line.strip()).union({0})

    # get reads to remove
    remove = set()
    with open(args.centrifuge) as f_centrifuge:
        try:
            header = next(f_centrifuge).strip().split('\t')
        except StopIteration:
            pass
        else:
            i_readid = header.index('readID')
            i_taxid = header.index('taxID')
            i_score = header.index('score')
            for line in f_centrifuge:
                try:
                    cols = line.strip().split('\t')
                    score = float(cols[i_score])
                    taxid = int(cols[i_taxid])
                    if taxid not in taxids_virus and score >= args.min_score:
                        remove.add(cols[i_readid])
                except:
                    # sometimes columns are not properly formatted (e.g. taxids with a .1)
                    # so we skip and simply keep those reads
                    continue

    # filter the reads
    with open(args.out, 'w') as out_f:
        for record in SeqIO.parse(args.fastq, 'fastq'):
            if record.id not in remove:
                SeqIO.write(record, out_f, 'fastq')


if __name__ == '__main__':
    main()
