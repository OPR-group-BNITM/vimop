#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.


import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='')
    parser.add_argument('--taxon-table', required=True, help='')
    parser.add_argument('--families', required=True, help='')
    parser.add_argument('--species', required=True, help='')
    args = parser.parse_args()

    first_seq = next(SeqIO.parse(args.input, 'fasta'))
    n_cols = len(str(first_seq.description).split('|'))

    if n_cols < 3:
        raise RuntimeError(
            "Expecting at least four columns in fasta header separated by '|' "
            "(Sequence-ID, description, family, species). Empty columns are allowed"
        )

    all_families = set()
    all_species = set()
    with open(args.taxon_table, 'w') as f_taxon_table:
        for i, seq in enumerate(SeqIO.parse(args.input, 'fasta'), 1):
            cols = str(seq.description).split('|')
            if len(cols) != n_cols:
                raise RuntimeError(f"Mismatch in number of columns for entry {i}")
            seqid, kingdom, family, species, *_ = cols
            f_taxon_table.write(f"{seqid}\t{kingdom}\t{family}\t{species}\n")
            all_families.add(family)
            all_species.add(species)
    with open(args.families, 'w') as f_fam:
        f_fam.write("\n".join(sorted(all_families)))
    with open(args.species, 'w') as f_species:
        f_species.write("\n".join(sorted(all_species)))


if __name__ == "__main__":
    main()
