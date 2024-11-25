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
        help="csv files with argets by filtering method (need to match assembly modes).",
        nargs='*'
    )
    parser.add_argument(
        "--virus-db-config",
        help="Config file of the virus data base (yaml)."
    )
    return parser


def read_unique_blast_hits(*fnames):
    """Get all targets found with blast but drop the duplicates.
    """
    cols = ['Reference', 'Description', 'Family', 'Organism', 'Segment', 'Orientation']
    blast_hits = pd.concat([
        pd.read_csv(fname, sep=',')[cols]
        for fname in fnames
    ]).drop_duplicates()
    return blast_hits


def rename_columns_extend_mapstats(df_mapping_stats):
    df_mapstats = df_mapping_stats.copy()
    df_mapstats.fillna({'ConsensusLength': 0, 'NCount': 0, 'Family': ''}, inplace=True)
    df_mapstats['CalledNucleobases'] = df_mapstats['ConsensusLength'] - df_mapstats['NCount']
    df_mapstats['Coverage'] = np.where(
        df_mapstats['ConsensusLength'] == 0,
        0,
        df_mapstats['CalledNucleobases'] / df_mapstats['ConsensusLength'] * 100
    ).round(0).astype(int)
    df_mapstats.rename(
        inplace=True,
        columns={
            'ReferenceLength': 'Length',
            'NCount': 'Ambiguous positions',
            'NumberOfMappedReads': 'Mapped reads',
            'AverageCoverage': 'Average read coverage',
            'CalledNucleobases': 'Positions called',
        }
    )
    return df_mapstats


def read_mapping_stats(fname_mapstats, fnames_blast_hits, virus_db_config):
    unique_blast_hits = read_unique_blast_hits(*fnames_blast_hits)
    mapping_stats = pd.read_csv(fname_mapstats, sep='\t')
    # Removing the '.1', '.2', etc. from the 'Reference' column
    mapping_stats['Reference'] = mapping_stats['Reference'].str.replace(r'\.\d+$', '', regex=True)
    mapping_stats = pd.merge(mapping_stats, unique_blast_hits, on='Reference')
    mapping_stats = rename_columns_extend_mapstats(mapping_stats)

    curated_organisms = {
        org.lower(): feats['name']
        for _, feats in virus_db_config['curated'].items()
        for org in feats['organisms']
    }
    mapping_stats['Curated'] = mapping_stats['Organism'].str.lower().isin(curated_organisms.keys())
    mapping_stats['Organism'] = mapping_stats['Organism'].map(curated_organisms).fillna(mapping_stats['Organism'])

    return mapping_stats


def html_report(
        df_mapping_stats,
        df_clean_read_stats,
        virus_db_config,
        fname_out 
):
    report = labs.BasicReport(report_title="Virus metagenommics sequencing")
    with report.add_section("Read Statistics", "Read Statistics"):
        df_readcounts = df_clean_read_stats.rename(columns={
            'step': 'Stage',
            'num_seqs': 'Reads passing the filter',
            'sum_len': 'Total Nucleotides',
            'min_len': 'Shortest Read',
            'avg_len': 'Mean Length',
            'max_len': 'Longest Read'
        })
        DataTable.from_pandas(
            df_readcounts[['Stage', 'Reads passing the filter', 'Mean Length']],
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
        feats['name'].lower(): feats['segments'] 
        for _, feats in virus_db_config['curated'].items()
    }
    for organism, mapstats in mapstats_curated.groupby("Organism"):
        section_name = organism
        with report.add_section(section_name, section_name):
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
                'Organism': organism,
                'Length': 0,
                'Coverage': 0,
                'Description': 'No hit found',
                'Positions called': 0,
                'Ambiguous positions': 0,
                'Mapped reads': 0,
                'Average read coverage': 0,
            }
            missing_segments = sorted(set(segments[organism.lower()]) - set(mapstats_segments['Segment'].values))
            missing_segments_rows = pd.DataFrame([{**{'Segment': seg}, **missing_segment_defaults} for seg in missing_segments])
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
                DataTable.from_pandas(mapstats_segments[cols_details], use_index=False, export=True)
                p(
                    """
                    Targets used for consensus building.
                    Coverage is reported in percent.
                    """
                )

    section_name = "Non-Curated"
    with report.add_section(section_name, section_name):
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
            DataTable.from_pandas(df_mapstats[cols], use_index=False, export=True)
            p(
                """
                Targets used for consensus building.
                Coverage is reported in percent.
                """
            )
    report.write(fname_out)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")

    with open(args.virus_db_config) as f_config:
        virus_db_config = yaml.safe_load(f_config)

    os.makedirs(args.outdir, exist_ok=True)
    df_mapping_stats = read_mapping_stats(args.mapping_stats, args.blast_hits, virus_db_config)
    df_mapping_stats.to_csv(os.path.join(args.outdir, f'{args.prefix}_consensus_stats.tsv'), sep='\t')
    df_clean_read_stats = pd.read_csv(args.clean_read_stats, sep='\t')

    html_report(
        df_mapping_stats,
        df_clean_read_stats,
        virus_db_config,
        os.path.join(args.outdir, f'{args.prefix}.html')
    )

    # TODO
    # - include assembly report
    # - per filter
    # -- number of reads
    # -- parameters chosen (min_readlen, min_)
    # -- succesful or not
    # -- number of contigs found
