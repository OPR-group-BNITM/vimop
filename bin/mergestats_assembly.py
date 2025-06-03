#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser('Merge assembly statisistics')
    parser.add_argument(
        '--stats',
        help='Input .tsv files.',
        required=True,
        nargs='*'
    )
    parser.add_argument(
        '--out',
        help='Output .tsv file.',
        required=True
    )
    args = parser.parse_args()

    assemblies = pd.concat([
        pd.read_csv(fname, sep='\t')[['step', 'num_seqs', 'avg_len']]
        for fname in args.stats
    ])

    split_cols = assemblies['step'].str.split('_', n=1, expand=True)
    assemblies['Sequence type'] = split_cols[0]
    assemblies['Stage'] = split_cols[1]

    assemblies = assemblies.rename(columns={
        'num_seqs': 'Reads/Contigs',
        'avg_len': 'Mean length',
    })

    assemblies.to_csv(args.out, sep="\t")


if __name__ == '__main__':
    main()
