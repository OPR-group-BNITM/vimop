from .util import get_named_logger, wf_parser  # noqa: ABS101
import pandas as pd
import yaml

from dominate.tags import h5, p, span, table, tbody, td, th, thead, tr
import ezcharts as ezc
from ezcharts.components.reports import labs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.snippets import Tabs


def argparser():
    parser = wf_parser("report_sample")
    parser.add_argument(
        '--virus-db-config',
        help='Config file of the virus data base (yaml).',
        required=True,
    )
    parser.add_argument(
        '--consensus-stats',
        help='.tsv file with consensus and mapping statistics.',
        required=True,
    )
    parser.add_argument(
        '--reads-stats',
        help='.tsv file with read statistics',
        required=True,
    )
    parser.add_argument(
        '--contigs-stats',
        help='.tsv file with contigs statistics',
        required=True,
    )
    parser.add_argument(
        '--out',
        help='Output filename',
        default='report.html'
    )
    return parser


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

    mapstats_curated = df_mapping_stats[df_mapping_stats['IsBest'] == True]
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
                mapstats_segments = mapstats[cols_all]
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


def read_tsv(fname):
    return pd.read_csv(fname, sep='\t')


def main(args):
    logger = get_named_logger("Report")

    # read input
    with open(args.virus_db_config) as f_config:
        virus_db_config = yaml.safe_load(f_config)

    consensus_stats = read_tsv(args.consensus_stats)
    clean_read_stats = read_tsv(args.reads_stats)
    contigs_stats = read_tsv(args.contigs_stats)

    # html report
    html_report(
        consensus_stats,
        clean_read_stats,
        contigs_stats,
        virus_db_config,
        args.out,
    )
