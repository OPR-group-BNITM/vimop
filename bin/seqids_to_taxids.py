#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.


import argparse


def read_taxid_table(fname_taxids):
    with open(fname_taxids) as f:
        return dict(zip(l.strip().split('\t') for l in f if l.strip()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxon-table', required=True, help='Family and species to each sequence ID.')
    parser.add_argument('--kingdom-taxids', required=True, help='Table with kingdom to taxID.')
    parser.add_argument('--species-taxids', required=True, help='Table with species to taxID.')
    parser.add_argument('--family-taxids', required=True, help='Table with family to taxID.')
    parser.add_argument('--seqid-to-taxid', required=True, help='Name of the output file.')
    args = parser.parse_args()

    taxids_species = read_taxid_table(args.species_taxids)
    taxids_families = read_taxid_table(args.family_taxids)
    taxids_kingdoms = read_taxid_table(args.kingdom_taxids)
    allowed_kingdom_str = ' '.join(taxids_kingdoms.keys())

    with open(args.taxon_table) as f_taxontable, open(args.seqid_to_taxid, 'w') as f_seqid_to_taxid:
        for line in enumerate(f_taxontable, 1):
            seqid, kingdom, family, species = line.strip().split('\t')
            if species in taxids_species:
                taxid = taxids_species[species]
            elif family in taxids_families:
                taxid = taxids_families[family]
            elif kingdom in taxids_kingdoms:
                taxid = taxids_kingdoms[kingdom]
            else:
                raise RuntimeError(
                    f"Sequence {seqid} does not have valid kingdom ({kingdom}), "
                    f"the following are allowed: {allowed_kingdom_str}"
                )
            f_seqid_to_taxid.write(f'{seqid}\t{taxid}\n')


if __name__ == "__main__":
    main()
