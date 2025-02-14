from .util import get_named_logger, wf_parser  # noqa: ABS101
import pandas as pd
import numpy as np


def argparser():
    parser = wf_parser('mergestats_read_clean')
    parser.add_argument(
        '--clean-read-stats',
        help='tsv file with read counts along the different cleaning steps.',
        required=True,
    )
    parser.add_argument(
        '--out',
        help='tsv file to write transformed data to.',
        default='out.tsv',
    )
    return parser


def read_clean_stats(fname):
    stats = pd.read_csv(fname, sep='\t')[['step', 'num_seqs', 'sum_len']]
    n_rows, _ = stats.shape
    entries = [
        {
            'Stage': 'Unfiltered',
            'Reads': stats.loc[0]['num_seqs'],
            'Nucleobases': stats.loc[0]['sum_len'],
        },
        *[
            {
                'Stage': stats.loc[i]['step'],
                'Reads': stats.loc[i-1]['num_seqs'] - stats.loc[i]['num_seqs'],
                'Nucleobases': stats.loc[i-1]['sum_len'] - stats.loc[i]['sum_len'],
            }
            for i in range(1, n_rows)
        ],
        {
            'Stage': 'Filtered',
            'Reads': stats.loc[n_rows-1]['num_seqs'],
            'Nucleobases': stats.loc[n_rows-1]['sum_len'],
        }
    ]
    stats_new = pd.DataFrame(entries)
    stats_new['Mean length'] = np.where(
        stats_new['Reads'] == 0,
        0,
        stats_new['Nucleobases'] / stats_new['Reads']
    )
    return stats_new


def main(args):
    logger = get_named_logger("Merge read statistics")
    clean_read_stats = read_clean_stats(args.clean_read_stats)
    clean_read_stats.to_csv(args.out, sep='\t')
