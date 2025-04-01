#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details


import argparse
import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml', required=True, help='yaml config file')
    parser.add_argument('--keys', required=True, help='Chain of keys', nargs='+')
    args = parser.parse_args()

    with open(args.yaml) as f_yaml:
        config = yaml.safe_load(f_yaml)

    retval = config
    for key in args.keys:
        retval = retval[key]
    
    print(retval)


if __name__ == '__main__':
    main()
