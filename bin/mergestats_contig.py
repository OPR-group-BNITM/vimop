#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
import pandas as pd


def df_read_and_merge(fnames, sep='\t'):
    return pd.concat([pd.read_csv(fname, sep=sep) for fname in fnames])


def merge_contigs_blasthits(contig_infos, blast_hits, contig_classes):

    cols_contigs = ['Contig', 'WorkflowMappingTarget', 'len', 'reads']
    contigs = contig_infos[cols_contigs]
    
    contigs= contigs.merge(contig_classes, left_on='Contig', right_on='readID', how='left')
    contigs = contigs.drop(columns=['readID'])

    # Merge the contig info with the blast hits
    merged = pd.merge(
        contigs,
        blast_hits,
        how='left',
        left_on=['Contig', 'WorkflowMappingTarget'],
        right_on=['Query', 'WorkflowMappingTarget'],
        suffixes=('_contig', '_blast')
    )

    # Compute coverages and sequence identity
    merged['Contig alignment coverage'] = ((merged['QueryFrom'] - merged['QueryTo']).abs() + 1) / merged['len']
    merged['Reference alignment coverage'] = ((merged['HitFrom'] - merged['HitTo']).abs() + 1) / merged['HitLength']
    merged['Sequence Identity'] = merged['IdenticalPositions'] / merged['AlignmentLength']

    # Rename and filter columns
    cols = {
        'WorkflowMappingTarget': 'Filter',
        'Contig': 'Contig',
        'name': 'Classification',
        'taxRank': 'Taxonomic Rank',
        'len': 'Length',
        'reads': 'Number of reads',
        'Reference': 'Blast Hit',
        'Organism': 'Organism',
        'HitLength': 'Hit length',
        'Contig alignment coverage': 'Contig alignment coverage',
        'Hit alignment coverage': 'Reference alignment coverage',
        'Sequence Identity': 'Sequence Identity',
    }
    merged.rename(columns=cols, inplace=True)
    # Fill empty values for contigs with no Blast hits
    merged.fillna(
        {   
            'Classification': 'Unclassified',
            'Blast Hit': 'Not found',
            'Organism': '',
            'Hit length': 0,
            'Contig alignment coverage': 0,
            'Reference alignment coverage': 0,
            'Sequence Identity': 0,
        },
        inplace=True
    )
    filtered = merged[cols.values()]
    filtered = filtered.drop_duplicates()

    return filtered.astype({'Hit length': int})


def main():

    parser = argparse.ArgumentParser('Merge contigs and blast info')
    parser.add_argument(
        '--contig-classes',
        help='csv files with centrifuge classification for each contigs by assembly mode.',
        nargs='*'
    )
    parser.add_argument(
        '--blast-hits',
        help='csv files with targets by filtering method (need to match assembly modes).',
        nargs='*'
    )
    parser.add_argument(
        '--contig-info',
        help='Tab-separated file with contig stats.',
        nargs='*'
    )
    parser.add_argument(
        '--out',
        help='.tsv file to write transformed data to.',
        default='out.tsv',
    )
    args = parser.parse_args()

    blast_hits = df_read_and_merge(args.blast_hits, sep=',')
    contigs = df_read_and_merge(args.contig_info)
    contig_classes = df_read_and_merge(args.contig_classes, sep=',')

    contigs_stats = merge_contigs_blasthits(contigs, blast_hits, contig_classes)
    
    contigs_stats.to_csv(args.out, sep='\t')

if __name__ == '__main__':
    main()
