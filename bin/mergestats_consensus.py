#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
import pandas as pd
import numpy as np
import yaml


def merge_mapstats_reference_info(mapstats, reference_info, virus_db_config):

    # Reduce blast hits to unique hits
    reference_info_cols = ['Reference', 'Description', 'Family', 'Organism', 'Segment', 'Orientation']
    unique_blast_hits = reference_info[reference_info_cols].drop_duplicates()

    merged = mapstats.copy()
    merged = pd.merge(merged, unique_blast_hits, on='Reference')

    # add information about curated virus database entries
    curated_organisms_labels = {
        org.lower(): organism_label
        for organism_label, feats in virus_db_config['curated'].items()
        for org in feats['organisms']
    }

    merged['Curated'] = merged['Organism'].str.lower().isin(curated_organisms_labels.keys())
    merged['Organism Label'] = merged['Organism'].str.lower().map(curated_organisms_labels).fillna('Non-Curated')

    # Fill missing values
    merged.fillna({'ConsensusLength': 0, 'NCount': 0, 'Family': ''}, inplace=True)
    
    # Add columns about called nucleobases and coverage
    merged['CalledNucleobases'] = merged['ConsensusLength'] - merged['NCount']
    merged['Coverage'] = np.where(
        merged['ConsensusLength'] == 0,
        0,
        merged['CalledNucleobases'] / merged['ConsensusLength'] * 100
    ).round(2)

    # Get the best representative for each curated organism + segment combination
    merged['IsBest'] = False
    max_indices = merged[merged['Curated']].groupby(['Organism Label', 'Segment'])['CalledNucleobases'].idxmax()
    merged.loc[max_indices, 'IsBest'] = True

    # rename columns
    merged.rename(
        inplace=True,
        columns={
            'ReferenceLength': 'Length',
            'NCount': 'Ambiguous positions',
            'NumberOfMappedReads': 'Mapped reads',
            'AverageCoverage': 'Average read coverage',
            'CalledNucleobases': 'Positions called',
        }
    )

    return merged


def main():

    parser = argparse.ArgumentParser("mergestats_consensus")
    parser.add_argument(
        '--virus-db-config',
        help='Config file of the virus data base (yaml).',
        required=True,
    )
    parser.add_argument(
        '--mapping-stats',
        help='tsv file with mapping statistics for each target.'
    )
    parser.add_argument(
        "--reference-info",
        help='.tsv file with information about the reference genomes.',
        required=True,
    )
    parser.add_argument(
        '--out',
        help='.tsv file to write transformed data to.',
        default='out.tsv',
    )
    args = parser.parse_args()

    with open(args.virus_db_config) as f_config:
        virus_db_config = yaml.safe_load(f_config)

    reference_info = pd.read_csv(args.reference_info, sep='\t').fillna('')
    mapstats = pd.read_csv(args.mapping_stats, sep='\t')

    consensus_stats = merge_mapstats_reference_info(
        mapstats, reference_info, virus_db_config   
    )
    consensus_stats.to_csv(args.out, sep='\t')


if __name__ == '__main__':
    main()
