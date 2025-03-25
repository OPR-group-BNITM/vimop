#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.

import argparse
import pandas as pd
import yaml
import os

from typing import List, Type
import base64

from bokeh.models import Title
from dominate.tags import (
    h1, p, div, section, html_tag,
    a, button, img, ul, li, strong,
    span
)

import ezcharts as ezc

from ezcharts.layout.snippets import Tabs, Grid
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.snippets.section import Section
from ezcharts.layout.snippets.tabs import ITabsClasses, ITabsStyles
from ezcharts.layout.snippets.banner import (
    IBannerStyles, IBannerClasses
)

from ezcharts.layout.base import Snippet
from ezcharts.layout.util import cls, css

from ezcharts.components.reports import Report, labs
from ezcharts.components.ezchart import EZChart
from ezcharts.components.theme import LAB_body_resources, LAB_head_resources


BACKGROUND_COLOR = '#1d1d1d'
PLOT_COLOR = None


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


class OprBanner(Snippet):
    """A styled div tag containing a heading and badges."""

    TAG = 'div'

    def __init__(self, report_title: str, pipeline_version: str) -> None:
        # TODO: replace the banner styles and classes?
        styles: IBannerStyles = IBannerStyles()
        classes: IBannerClasses = IBannerClasses()
        # Override the default background color
        styles.container = css(
            f"background-color: {BACKGROUND_COLOR} !important;",
            "padding: 15px;",
            "border-radius: 5px;"
        )

        super().__init__(
            styles=styles,
            classes=classes,
            className=classes.container,
            style=styles.container
        )
        with self:
            with div(className=self.classes.inner, style=self.styles.inner):
                h1(report_title)
                p(f'Pipeline version {pipeline_version}')
                p('Research use only')


class OprLogo(div):
    """OPR logo element."""
    def __init__(self) -> None:
        fname_logo = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'OPR_logo_v01-light.cropped.png'
        )

        # Convert image to Base64
        with open(fname_logo, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode("utf-8")

        super().__init__(
            img(src=f"data:image/png;base64,{encoded_string}", style="height: 75px;", alt="OPR Logo"),
            tagname='div',
            className="d-flex"
        )


class OprNavigation(Snippet):
    """A styled nav component for use in a Report."""

    TAG = 'nav'

    def __init__(
        self,
        logo: Type[html_tag],
        groups: List[str],
        header_height: int = 75,
        classes: labs.ILabsNavigationClasses = labs.ILabsNavigationClasses()
    ) -> None:
        spacer = div(
            className=classes.spacer,
            style=f"margin-top: {header_height}px;"
        )
        super().__init__(
            styles=None,
            classes=classes,
            style=f"min-height: {header_height}px; background-color: {BACKGROUND_COLOR} !important;",
            className=classes.container
        )
        spacer.add(self)
        with self:
            with div(className=self.classes.inner):
                with a(href="https://github.com/OPR-group-BNITM",
                       className=self.classes.logo):
                    logo()
                button(
                    "Jump to section... ",
                    cls=self.classes.dropdown_btn,
                    type="button",
                    id="dropdownMenuButton",
                    data_bs_toggle="dropdown",
                    aria_haspopup="true",
                    aria_expanded="false"
                )
                ngroups = len(groups)
                with div(className=self.classes.dropdown_menu):
                    for count, group in enumerate(groups):
                        setattr(
                            self, group,
                            div(className='', __pretty=False)
                        )
                        if count != ngroups - 1:
                            div(cls="dropdown-divider")

    def add_link(
        self,
        group: str,
        link_title: str,
        link_href: str
    ) -> None:
        """Add a header nav link to the header links list."""
        group_list = getattr(self, group)
        with group_list:
            a(
                link_title,
                href=link_href,
                className=self.classes.dropdown_item_link
            )


class OprReport(Report):
    """A basic OPR-themed report."""

    def __init__(self, title, pipeline_version):
        super().__init__(
            report_title=title,
            head_resources=LAB_head_resources,
            body_resources=LAB_body_resources
        )
        with self.header:
            self.header.attributes["style"] = f"background-color: {BACKGROUND_COLOR} !important;"
            self.nav = OprNavigation(logo=OprLogo, groups=['main', 'meta'])
            self.intro_content = section(
                id="intro-content",
                role="region",
                style=f"background-color: {BACKGROUND_COLOR}  !important;"
            )
            with self.intro_content:
                self.banner = OprBanner(title, pipeline_version)

        with self.main:
            self.main_content = section(id="main-content", role="region")

    def add_section(
        self,
        title: str,
        link: str,
        overflow: bool = False
    ) -> Section:
        """Add a section to the main_content region."""
        href = self.get_uid('Section')
        self.nav.add_link('main', link, f'#{href}')
        with self.main_content:
            return Section(href, title, overflow=overflow)


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
        color=PLOT_COLOR
    )
    subtitle = f"Mean: {mean:.1f}. Median: {median:.1f}"

    plot._fig.xaxis.axis_label = xaxis_label
    plot._fig.add_layout(Title(text=subtitle, text_font_size="0.8em"), 'above')
    plot._fig.add_layout(Title(text=title, text_font_size="1.5em"), 'above')
    
    plot._fig.x_range.start = min_val
    plot._fig.x_range.end = max_val

    plot._fig.yaxis.axis_label = 'Number of reads'

    return plot


class Legends:

    def __init__(self, descriptions: dict, default_title: str = None):
        self.descriptions = descriptions
        self.default_title = default_title if default_title else ""

    def legend(self, columns: list, title: str = None):
        with div(style="font-size: 0.9em; margin-top: 1em; color: #5c5c5c;"):
            title = title if title else self.default_title
            p(title, style="font-size: 1em; font-weight: bold; margin-bottom: 0.3em;")
            with ul(style="padding-left: 1.2em; margin: 0;"):
                for col in columns:
                    with li():
                        strong(f"{col}: ")
                        span(self.descriptions[col])


def html_report(
        pipeline_version,
        samplename,
        df_mapping_stats,
        df_clean_read_stats,
        df_contig_stats,
        df_assembly_stats,
        df_lenqual_trim,
        df_lenqual_clean,
        virus_db_config,
        fname_out 
):
    report = OprReport(
        title=f"Virus metagenomics sequencing report for sample {samplename}",
        pipeline_version=pipeline_version
    )
    with report.add_section("Read Statistics", "Read Statistics"):
        tabs_readstats = NoMarginTabs()
        with tabs_readstats.add_tab("Table"):
            cols = ['Stage', 'Reads', 'Mean length']
            DataTable.from_pandas(
                df_clean_read_stats[cols].round({'Mean length': 0}).astype({'Mean length': int}),
                use_index=False,
                export=False
            )
            Legends(
                {
                    'Stage': 'Host/contaminant filtering step',
                    'Reads': 'Number of reads in this group',
                    'Mean length': 'Mean read length within this read group',
                },
                'Reads in groups due to host/contaminant filtering'
            ).legend(cols)
        with tabs_readstats.add_tab("Distributions trimmed"):
            try:
                with Grid(columns=2):
                    EZChart(histogram_plot(df_lenqual_trim['Length'], 'Length', 'Base pairs'))
                    EZChart(histogram_plot(df_lenqual_trim['Quality'], 'Quality', 'Avgerage quality score per read'))
                p(
                    """
                    The distributions of the reads after trimming. They still contain host reads.
                    """
                )
            except ValueError:
                nseqs = len(df_lenqual_trim['Length'])
                p(f'Failed to create histograms for {nseqs} sequences')
        with tabs_readstats.add_tab("Distributions cleaned"):
            try:
                with Grid(columns=2):
                    EZChart(histogram_plot(df_lenqual_clean['Length'], 'Length', 'Base pairs'))
                    EZChart(histogram_plot(df_lenqual_clean['Quality'], 'Quality', 'Avgerage quality score per read'))
                p(
                    """
                    The distributions of the reads after cleaning. Host reads and technical contaminants are removed.
                    """
                )
            except ValueError:
                nseqs = len(df_lenqual_clean['Length'])
                p(f'Failed to create histograms for {nseqs} sequences')
    mapstats_curated = df_mapping_stats[df_mapping_stats['IsBest'] == True]
    segments = {
        label: feats['segments'] 
        for label, feats in virus_db_config['curated'].items()
    }

    section_name = 'Consensus'
    consensus_legends = Legends(
        {
            'Reference': 'ID of the reference sequence',
            'Family': 'Virus family',
            'Organism': 'Name of the species',
            'Segment': 'Identifier of the segment. Unsegmented for single segment virus genomes',
            'Length': 'Length of the reference genome',
            'Coverage': 'Percent of the reference genome that were succesfully called',
            'Description': 'Data base description of the reference',
            'Positions called': 'Number of bases called',
            'Ambiguous positions': 'Number of ambiguous positions set to "N"',
            'Mapped reads': 'Number of reads aligned to the reference genome',
            'Average read coverage': 'Average number of reads per reference genome position',
        },
        "Reference-based genome assembly statistics",
    )
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
                    'Description',
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
                    'Average read coverage',
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
                    consensus_legends.legend(cols_overview)
                with tabs.add_tab("Details"):
                    DataTable.from_pandas(mapstats_segments[cols_details].round(round), use_index=False, export=True)
                    consensus_legends.legend(cols_details)
        with tabs_consensus.add_tab("Non-Curated"):
            df_mapstats = df_mapping_stats[df_mapping_stats['Curated'] == False]
            tabs = NoMarginTabs()
            with tabs.add_tab("Overview"):
                cols = ['Reference', 'Family', 'Organism', 'Length', 'Coverage', 'Description']
                DataTable.from_pandas(df_mapstats[cols], use_index=False, export=True)
                consensus_legends.legend(cols)
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
                consensus_legends.legend(cols)
        with tabs_consensus.add_tab("All"):
            df_mapstats = df_mapping_stats
            tabs = NoMarginTabs()
            with tabs.add_tab("Overview"):
                cols = ['Reference', 'Family', 'Organism', 'Length', 'Coverage', 'Description']
                DataTable.from_pandas(df_mapstats[cols], use_index=False, export=True)
                consensus_legends.legend(cols)
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
                consensus_legends.legend(cols)

    section_name = "Contigs"
    with report.add_section(section_name, section_name):
        tabs_contigs = NoMarginTabs()
        with tabs_contigs.add_tab("Contigs"):
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
            Legends(
                {
                    'Filter': 'Filter used on reads before assembly',
                    'Contig': 'Contig identifier',
                    'Length': 'Length of the contig in base pairs',
                    'Number of reads': 'Number of (corrected) reads used by canu to build this contig',
                    'Blast Hit': 'Virus reference genome found with blast search',
                    'Organism': 'Organism of the blast hit',
                    'Hit length': 'Length of the blast hit reference genome in base pairs',
                    'Contig alignment coverage': 'Share of the contig aligned by blast',
                    'Reference alignment coverage': 'Share of the reference aligned by blast',
                    'Sequence Identity': 'Sequence similarity of the aligned parts',
                },
                'Contigs and targets found'
            ).legend(cols)
        with tabs_contigs.add_tab("Assembly Statistics"):
            cols = ['Stage', 'Sequence type', 'Reads/Contigs', 'Mean length']
            digits = {
                'Reads/Contigs': 0,
                'Mean length': 1,
            }
            DataTable.from_pandas(df_assembly_stats[cols], use_index=False, export=True)
            Legends(
                {
                    'Stage': 'Assembly run (e.g. for a given filter, no-filter or re-assemblies)',
                    'Sequence type': 'Type of the sequence',
                    'Reads/Contigs': 'Number of reads or contigs',
                    'Mean length': 'Mean length of reads or contigs in base pairs',
                },
                'Read and contig numbers in the different assembly runs'
            ).legend(cols)

    report.write(fname_out)


def read_tsv(fname):
    return pd.read_csv(fname, sep='\t').fillna('')


def main():
    parser = argparse.ArgumentParser("report_sample")
    parser.add_argument(
        '--pipeline-version',
        help='The version of this pipeline.',
        required=True,
    )
    parser.add_argument(
        '--samplename',
        help='The name of your sample.',
        required=True,
    )
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
        '--assembly-stats',
        help='.tsv file with assembly statistics',
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
    args = parser.parse_args()

    with open(args.virus_db_config) as f_config:
        virus_db_config = yaml.safe_load(f_config)

    consensus_stats = read_tsv(args.consensus_stats)
    clean_read_stats = read_tsv(args.reads_stats)
    contigs_stats = read_tsv(args.contigs_stats)
    assembly_stats = read_tsv(args.assembly_stats)

    lenqual_trim = read_tsv(args.trimmed_read_distribution)
    lenqual_clean = read_tsv(args.cleaned_read_distribution)

    html_report(
        args.pipeline_version,
        args.samplename,
        consensus_stats,
        clean_read_stats,
        contigs_stats,
        assembly_stats,
        lenqual_trim,
        lenqual_clean,
        virus_db_config,
        args.out,
    )


if __name__ == '__main__':
    main()
