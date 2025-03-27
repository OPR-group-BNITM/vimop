#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
import yaml
import pandas as pd


def read_yaml(fname):
    with open(fname) as f_config:
        return yaml.safe_load(f_config)


def main():
    
    parser = argparse.ArgumentParser("mergestats_consensus")
    parser.add_argument(
        '--virus-db-config',
        help='Config file of the virus data base (yaml).',
        required=True,
    )
    parser.add_argument(
        '--contamination-db-config',
        help='Config file of the contamination data base (yaml).',
        required=True,
    )
    parser.add_argument(
        '--classification-db-config',
        help='Config file of the classification data base (yaml).',
        required=True,
    )
    parser.add_argument(
        '--out',
        help='.tsv file to write transformed data to.',
        default='out.tsv',
    )
    args = parser.parse_args()

    configs = {
        'Filters': read_yaml(args.contamination_db_config),
        'Classification': read_yaml(args.classification_db_config),
        'Virus': read_yaml(args.virus_db_config),
    }
    data = [
        {
            'Data Base': db_name,
            'Version': config['version'],
            'Description': config['description'],

        }
        for db_name, config in configs.items()
    ]
    pd.DataFrame(data).to_csv(args.out, sep='\t')


if __name__ == '__main__':
    main()
