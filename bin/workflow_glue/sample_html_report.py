from .util import get_named_logger, wf_parser  # noqa: ABS101
import pandas as pd
import yaml

from bokeh.models import Title
from dominate.tags import h5, p, span, table, tbody, td, th, thead, tr
import ezcharts as ezc
from ezcharts.components.reports import labs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.snippets import Tabs, Grid
from ezcharts.components.ezchart import EZChart

# ezcharts wrappers
from ezcharts.layout.util import cls, css
from ezcharts.layout.snippets.tabs import ITabsClasses, ITabsStyles


class NoMarginTabsClasses(ITabsClasses):
    """Override tab classes to remove margins."""
    tab_buttons_list: str = cls("nav", "nav-tabs")  # Removed "mb-2"
    tab_contents_container: str = cls("tab-content", "p-0")  # Set padding to 0


class NoMarginTabsStyles(ITabsStyles):
    """Override tab styles to remove margins."""
    tab_button: str = css(
        "margin-bottom: 0",  # Remove negative margin
        "font-weight: 600",
        "cursor: pointer",
        "border-color: transparent"
    )
    tab_button_active: str = css(
        "border-bottom: 2px solid #0079a4",
        "color: #0079a4!important"
    )


class NoMarginTabs(Tabs):
    """Custom Tabs class without extra margins."""
    def __init__(self) -> None:
        super().__init__(styles=NoMarginTabsStyles(), classes=NoMarginTabsClasses())


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
        '--trimmed-read-distribution',
        help='.tsv file with lengths and qualities of the trimmed reads.',
        required=True,
    )
    parser.add_argument(
        '--cleaned-read-distribution',
        help='.tsv file with lengths and qualities of the trimmed reads.',
        required=True,
    )
    parser.add_argument(
        '--out',
        help='Output filename',
        default='report.html'
    )
    return parser


def histogram_plot(data, title, xaxis_label):

    nbins = 300
    min_val = data.min()
    max_val = data.max()
    bin_width = (max_val - min_val) / nbins
    mean = data.mean()
    median = data.median()

    plot = ezc.histplot(
        data=data,
        binwidth=bin_width,
        binrange=(min_val, max_val),
        color=None
    )
    subtitle = f"Mean: {mean:.1f}. Median: {median:.1f}"

    plot._fig.xaxis.axis_label = xaxis_label
    plot._fig.add_layout(Title(text=subtitle, text_font_size="0.8em"), 'above')
    plot._fig.add_layout(Title(text=title, text_font_size="1.5em"), 'above')
    
    plot._fig.x_range.start = min_val
    plot._fig.x_range.end = max_val

    plot._fig.yaxis.axis_label = 'Number of reads'

    return plot


def html_report(
        df_mapping_stats,
        df_clean_read_stats,
        df_contig_stats,
        df_lenqual_trim,
        df_lenqual_clean,
        virus_db_config,
        fname_out 
):
    report = labs.BasicReport(report_title="Virus metagenommics sequencing")
    with report.add_section("Read Statistics", "Read Statistics"):
        tabs_readstats = NoMarginTabs()
        with tabs_readstats.add_tab("Table"):
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
        with tabs_readstats.add_tab("Distributions trimmed"):
            with Grid(columns=2):
                EZChart(histogram_plot(df_lenqual_trim['Length'], 'Length', 'Base pairs'))
                EZChart(histogram_plot(df_lenqual_trim['Quality'], 'Quality', 'Avgerage quality score per read'))
            p(
                """
                The distributions of the reads after trimming. They still contain host reads.
                """
            )
        with tabs_readstats.add_tab("Distributions cleaned"):
            with Grid(columns=2):
                EZChart(histogram_plot(df_lenqual_clean['Length'], 'Length', 'Base pairs'))
                EZChart(histogram_plot(df_lenqual_clean['Quality'], 'Quality', 'Avgerage quality score per read'))
            p(
                """
                The distributions of the reads after cleaning. Host reads and technical contaminants are removed.
                """
            )
    mapstats_curated = df_mapping_stats[df_mapping_stats['IsBest'] == True]
    segments = {
        label: feats['segments'] 
        for label, feats in virus_db_config['curated'].items()
    }

    section_name = 'Consensus'
    with report.add_section(section_name, section_name):
        tabs_consensus = NoMarginTabs()
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
                tabs = NoMarginTabs()
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
            tabs = NoMarginTabs()
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
    return pd.read_csv(fname, sep='\t').fillna('')


def main(args):
    logger = get_named_logger("Report")

    with open(args.virus_db_config) as f_config:
        virus_db_config = yaml.safe_load(f_config)

    consensus_stats = read_tsv(args.consensus_stats)
    clean_read_stats = read_tsv(args.reads_stats)
    contigs_stats = read_tsv(args.contigs_stats)

    lenqual_trim = read_tsv(args.trimmed_read_distribution)
    lenqual_clean = read_tsv(args.cleaned_read_distribution)

    html_report(
        consensus_stats,
        clean_read_stats,
        contigs_stats,
        lenqual_trim,
        lenqual_clean,
        virus_db_config,
        args.out,
    )
