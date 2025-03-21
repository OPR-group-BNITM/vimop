import argparse
import os
import shutil
import pandas as pd


def reference_id_from_fasta(fname):
    with open(fname) as f:
        return next(f).split('reference=')[1].split()[0]


def main():

    parser = argparse.ArgumentParser("mergestats_consensus")
    parser.add_argument(
        '--consensus-stats',
        help='tsv file with consensus statistics.',
        required=True,
    )
    parser.add_argument(
        '--consensus-files',
        help='fasta files with consensus sequences.',
        nargs='*',
    )
    parser.add_argument(
        '--out-dir',
        help='.tsv file to write transformed data to.',
        default='out.tsv',
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    seqs_by_id = {
        reference_id_from_fasta(fname): fname
        for fname in args.consensus_files
    }
    stats = pd.read_csv(args.consensus_stats, sep='\t')
    best = stats[stats['IsBest'] == True]

    for _, row in best.iterrows():

        refid = row['Reference']
        label = row['Organism Label']
        segment = row['Segment']

        fname_out = (
            f'{label}.fasta'
            if segment == 'Unsegmented'
            else f'{label}_{segment}.fasta'
        )
        shutil.copyfile(
            seqs_by_id[refid],
            os.path.join(args.out_dir, fname_out)
        )


if __name__ == '__main__':
    main()
