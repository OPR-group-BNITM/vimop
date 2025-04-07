#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.


# checksum_tool_zipped: sha256sum
# checksum_tool_directories: "find . -type f -exec sha256sum {} \; | sort | sha256sum | cut -d ' ' -f1"


import argparse
import os
import subprocess
import urllib.request
import yaml


def read_yaml(fname):
    with open(fname) as f_config:
        return yaml.safe_load(f_config)


def compute_directory_checksum(directory):
    result = subprocess.run(
        "find . -type f -exec sha256sum {} \\; | sort | sha256sum | cut -d ' ' -f1",
        shell=True,
        cwd=directory,
        stdout=subprocess.PIPE,
        text=True
    )
    return result.stdout.strip()


def file_checksum(fname):
    result = subprocess.run(
        f"sha256sum {fname} | cut -d ' ' -f1",
        shell=True,
        stdout=subprocess.PIPE,
        text=True
    )
    return result.stdout.strip()


def check_sum(directory, checksum_to_check):
    return checksum_to_check == compute_directory_checksum(directory)


def check_update_required(directory, checksum_to_check):
    if not os.path.isdir(directory):
        return True
    return not check_sum(directory, checksum_to_check)


def download_db(url, fname_out, checksum_tar):
    urllib.request.urlretrieve(url, fname_out)
    checksum_dl = file_checksum(fname_out)
    if checksum_tar != checksum_dl:
        raise RuntimeError(
            f"sha256 sum of {fname_out} is {checksum_dl} does not "
            f"match configuration ({checksum_tar})"
        )


def extract_tar_xz(fname, output_dir="."):
    subprocess.run(
        ["tar", "-xf", fname, "-C", output_dir],
        check=True
    )


def extract_db(fname, expected_checksum):
    extract_tar_xz(fname)
    outdir_name = fname.split('.')[0]
    directory_checksum = compute_directory_checksum(outdir_name)
    if directory_checksum == expected_checksum:
        raise RuntimeError(
            f'Checksum of extracted directory ({directory_checksum}) '
            f'does not match expected value ({expected_checksum})'
        )


def main_OLD():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--db-basedir',
        help='',
        required=True,
    )
    parser.add_argument(
        '--db-update-config',
        help='',
        required=True,
    )
    parser.add_argument(
        '--to-update',
        help='',
        required=True,
    )
    args = parser.parse_args()

    update_config = read_yaml(args.db_update_config)
    possible_cases = update_config['sub_databases']

    to_update = possible_cases if args.to_update.lower() == 'all' else args.to_update.split(',')

    for case in to_update:
        if case not in possible_cases:
            raise RuntimeError(f'{case} not in {",".join(sorted(possible_cases))}')

    update_required = [
        case
        for case in to_update
        if check_update_required(
            os.path.join(args.db_basedir, case),
            update_config['sub_databases'][case]['checksum_directory']
        )
    ]

    for case in update_required:
        conf = update_config['sub_databases'][case]
        fname_tarxz = f'{case}.tar.xz'
        download_db(conf['url'], fname_tarxz, conf['checksum_zipped'])
        extract_db(fname_tarxz, conf['checksum_directory'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--db-basedir',
        help='',
        required=True,
    )
    parser.add_argument(
        '--db-update-config',
        help='',
        required=True,
    )
    parser.add_argument(
        '--to-update',
        help='',
        required=True,
    )
    args = parser.parse_args()
    


if __name__ == '__main__':
    main()
