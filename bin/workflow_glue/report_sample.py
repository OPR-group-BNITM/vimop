from .util import get_named_logger, wf_parser  # noqa: ABS101
import pandas as pd
import os
import numpy as np
import yaml

from dominate.tags import h5, p, span, table, tbody, td, th, thead, tr
import ezcharts as ezc
from ezcharts.components.reports import labs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.snippets import Tabs


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report_sample")
    parser.add_argument("--prefix", help="Filename prefix.", default='sample_report')
    parser.add_argument("--outdir", help="Output directory.", default='.')
    parser.add_argument(
        "--mapping-stats",
        help="tsv file with mapping statistics for each target."
    )
    parser.add_argument(
        "--clean-read-stats",
        help="tsv file with read counts along the different cleaning steps."
    )
    parser.add_argument(
        '--assembly-modes',
        help='Names of the different assembly modes.',
        nargs='*'
    )
    parser.add_argument(
        "--assembly-read-stats",
        help="tsv files with read counts along the different assembly runs (need to match assembly modes).",
        nargs='*'
    )
    parser.add_argument(
        "--blast-hits",
        help="csv files with targets by filtering method (need to match assembly modes).",
        nargs='*'
    )
    parser.add_argument(
        "--contig-info",
        help="Tab-separated file with contig stats.",
        nargs='*'
    )
    parser.add_argument(
        "--virus-db-config",
        help="Config file of the virus data base (yaml)."
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


def df_read_and_merge(fnames, sep='\t'):
    return pd.concat([pd.read_csv(fname, sep=sep) for fname in fnames])


def merge_mapstats_blasthits(mapstats, blasthits, virus_db_config):

    # Reduce blast hits to unique hits
    blast_cols = ['Reference', 'Description', 'Family', 'Organism', 'Segment', 'Orientation']
    unique_blast_hits = blasthits[blast_cols].drop_duplicates()

    merged = mapstats.copy()

    # Removing the '.1', '.2', etc. from the 'Reference' column and then merge
    merged['Reference'] = merged['Reference'].str.replace(r'\.\d+$', '', regex=True)
    merged = pd.merge(merged, unique_blast_hits, on='Reference')

    # add information about curated virus database entries
    curated_organisms_labels = {
        org.lower(): organism_label
        for organism_label, feats in virus_db_config['curated'].items()
        for org in feats['organisms']
    }
    curated_organisms_names = {
        label: feats['name']
        for label, feats in virus_db_config['curated'].items()
    }

    merged['Curated'] = merged['Organism'].str.lower().isin(curated_organisms_labels.keys())
    merged['Organism Label'] = merged['Organism'].str.lower().map(curated_organisms_labels).fillna('Non-Curated')
    merged['Organism'] = merged['Organism'].map(curated_organisms_names).fillna(merged['Organism'])

    # Fill missing values
    merged.fillna({'ConsensusLength': 0, 'NCount': 0, 'Family': ''}, inplace=True)
    
    # Add columns about called nucleobases and coverage
    merged['CalledNucleobases'] = merged['ConsensusLength'] - merged['NCount']
    merged['Coverage'] = np.where(
        merged['ConsensusLength'] == 0,
        0,
        merged['CalledNucleobases'] / merged['ConsensusLength'] * 100
    ).round(0).astype(int)

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


def merge_contigs_blasthits(contig_infos, blast_hits):

    cols_contigs = ['Contig', 'WorkflowMappingTarget', 'len', 'reads']
    contigs = contig_infos[cols_contigs]

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

    return filtered.astype({'Hit length': int})


def html_report(
        df_mapping_stats,
        df_clean_read_stats,
        df_contig_stats,
        virus_db_config,
        fname_out 
):
    report = labs.BasicReport(report_title="Virus metagenommics sequencing")
    with report.add_section("Read Statistics", "Read Statistics"):
        DataTable.from_pandas(
            df_clean_read_stats[['Stage', 'Reads', 'Mean length']].round({'Mean length': 0}).astype({'Mean length': int}),
            use_index=False,
            export=False
        )
        p(
            """
            Reads left after each filtering step.
            """
        )

    mapstats_curated = df_mapping_stats[df_mapping_stats['Curated'] == True]
    segments = {
        label: feats['segments'] 
        for label, feats in virus_db_config['curated'].items()
    }

    section_name = 'Consensus'
    with report.add_section(section_name, section_name):
        tabs_consensus = Tabs()
        for organism_label, mapstats in mapstats_curated.groupby("Organism Label"):
            with tabs_consensus.add_tab(organism_label):
                cols_overview = [
                    'Reference',
                    'Family',
                    'Organism',
                    'Segment',
                    'Length',
                    'Coverage',
                    'Description'
                ]
                cols_details = [
                    'Reference',
                    'Family',
                    'Organism',
                    'Segment',
                    'Length',
                    'Coverage',
                    'Positions called',
                    'Ambiguous positions',
                    'Mapped reads',
                    'Average read coverage'
                ]
                cols_all = list(set(cols_overview + cols_details))
                mapstats_segments = mapstats.loc[mapstats.groupby('Segment')['Positions called'].idxmax()][cols_all]
                mapstats_segments.reset_index(drop=True, inplace=True)
                missing_segment_defaults = {
                    'Reference': 'Not found',
                    'Family': '',
                    'Organism': virus_db_config['curated'][organism_label]['name'],
                    'Length': 0,
                    'Coverage': 0,
                    'Description': 'No hit found',
                    'Positions called': 0,
                    'Ambiguous positions': 0,
                    'Mapped reads': 0,
                    'Average read coverage': 0,
                }
                round = {'Average read coverage': 1}
                missing_segments = sorted(set(segments[organism_label]) - set(mapstats_segments['Segment'].values))
                missing_segments_rows = pd.DataFrame([
                    {**{'Segment': seg}, **missing_segment_defaults}
                    for seg in missing_segments
                ])
                mapstats_segments = pd.concat(
                    [mapstats_segments, missing_segments_rows],
                    ignore_index=True
                )
                tabs = Tabs()
                with tabs.add_tab("Overview"):
                    DataTable.from_pandas(mapstats_segments[cols_overview], use_index=False, export=True)
                    p(
                        """
                        Targets used for consensus building.
                        Coverage is reported in percent.
                        """
                    )
                with tabs.add_tab("Details"):
                    DataTable.from_pandas(mapstats_segments[cols_details].round(round), use_index=False, export=True)
                    p(
                        """
                        Targets used for consensus building.
                        Coverage is reported in percent.
                        """
                    )
        with tabs_consensus.add_tab("Non-Curated"):
            df_mapstats = df_mapping_stats[df_mapping_stats['Curated'] == False]
            tabs = Tabs()
            with tabs.add_tab("Overview"):
                cols = ['Reference', 'Family', 'Organism', 'Length', 'Coverage', 'Description']
                DataTable.from_pandas(df_mapstats[cols], use_index=False, export=True)
                p(
                    """
                    Targets used for consensus building.
                    Coverage is reported in percent.
                    """
                )
            with tabs.add_tab("Details"):
                cols = [
                    'Reference',
                    'Family',
                    'Organism',
                    'Length',
                    'Coverage',
                    'Positions called',
                    'Ambiguous positions',
                    'Mapped reads',
                    'Average read coverage'
                ]
                round = {'Average read coverage': 1}
                DataTable.from_pandas(df_mapstats[cols].round(round), use_index=False, export=True)
                p(
                    """
                    Targets used for consensus building.
                    Coverage is reported in percent.
                    """
                )

    section_name = "Contigs"
    with report.add_section(section_name, section_name):
        # TODO: add tabs for overview and table?
        cols = [
            'Filter',
            'Contig',
            'Length',
            'Number of reads',
            'Blast Hit',
            'Organism',
            'Hit length',
            'Contig alignment coverage',
            'Reference alignment coverage',
            'Sequence Identity',
        ]
        digits = {
            'Contig alignment coverage': 2,
            'Reference alignment coverage': 2,
            'Sequence Identity': 2,
        }
        DataTable.from_pandas(df_contig_stats.round(digits)[cols], use_index=False, export=True)
    report.write(fname_out)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")

    os.makedirs(args.outdir, exist_ok=True)

    # read input
    with open(args.virus_db_config) as f_config:
        virus_db_config = yaml.safe_load(f_config)
    blast_hits = df_read_and_merge(args.blast_hits, sep=',')
    contigs = df_read_and_merge(args.contig_info)
    mapstats = pd.read_csv(args.mapping_stats, sep='\t')
    clean_read_stats = read_clean_stats(args.clean_read_stats)

    # merge tables
    consensus_stats = merge_mapstats_blasthits(mapstats, blast_hits, virus_db_config)
    contigs_stats = merge_contigs_blasthits(contigs, blast_hits)

    # write output
    consensus_stats.to_csv(
        os.path.join(args.outdir, f'{args.prefix}_consensus_stats.tsv'),
        sep='\t'
    )

    # html report
    fname_html_out = os.path.join(args.outdir, f'{args.prefix}.html')
    html_report(
        consensus_stats,
        clean_read_stats,
        contigs_stats,
        virus_db_config,
        fname_html_out,
    )

    # TODO
    # - include assembly overview report
    # - per filter
    # -- number of reads input/used
    # -- parameters chosen (min_readlen, sampling?)
    # -- succesful or not
    # -- number of contigs found
