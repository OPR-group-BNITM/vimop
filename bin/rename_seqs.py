#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser("rename_seqs")
    parser.add_argument(
        '--prefix',
        help='Prefix to add before numbering',
        required=True,
    )
    parser.add_argument(
        '--input',
        help='Input fasta file.',
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Output fasta file',
        required=True,
    )
    args = parser.parse_args()

    with open(args.output, 'w') as f_out:
        for i, seq in enumerate(SeqIO.parse(args.input, 'fasta'), 1):
            try:
                nreads = int(seq.description.split("reads=")[1].split()[0])
            except:
                nreads = 1
            new_id = f'{args.prefix}{i}'
            description = f'len={len(seq.seq)} reads={nreads} original_name={seq.id}'
            new_record = SeqRecord(
                id=new_id,
                name=new_id,
                description=description,
                seq=seq.seq
            )
            SeqIO.write(new_record, f_out, 'fasta')


if __name__ == '__main__':
    main()
