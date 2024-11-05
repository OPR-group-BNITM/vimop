from .util import get_named_logger, wf_parser  # noqa: ABS101
import pandas as pd
import os

from dominate.tags import h5, p, span, table, tbody, td, th, thead, tr
import ezcharts as ezc
from ezcharts.components.reports import labs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.snippets import Tabs

# from pyecharts.charts import Page
# from pyecharts.components import Table
# from pyecharts.options import ComponentTitleOpts


# def html_table(df: pd.DataFrame, title: str) -> Table:
#     table = Table()
#     headers = df.columns.tolist()
#     rows = df.values.tolist()
#     table.add(headers, rows)
#     table.set_global_opts(title_opts=ComponentTitleOpts(title=title))
#     return table


# def html_report(
#         df_mapping_stats,
#         fname_out
# ):
#     page = Page()
#     page.add(html_table(df_mapping_stats, "Consensus stats"))
#     page.render(fname_out)


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
    return parser


def read_unique_blast_hits(*fnames):
    """Get all targets found with blast but drop the duplicates.
    """
    cols = ['Reference', 'Description', 'Family', 'Organism']
    blast_hits = pd.concat([
        pd.read_csv(fname, sep=',')[cols]
        for fname in fnames
    ]).drop_duplicates()
    return blast_hits


def read_mapping_stats(fname, fnames_blast_hits):
    unique_blast_hits = read_unique_blast_hits(*fnames_blast_hits)
    mapping_stats = pd.read_csv(fname, sep='\t')
    # Removing the '.1', '.2', etc. from the 'Reference' column
    mapping_stats['Reference'] = mapping_stats['Reference'].str.replace(r'\.\d+$', '', regex=True)
    mapping_stats_with_target_desciptions = pd.merge(mapping_stats, unique_blast_hits, on="Reference")
    return mapping_stats_with_target_desciptions


def html_report(
        df_mapping_stats,
        df_clean_read_stats,
        fname_out 
):
    report = labs.BasicReport(report_title="Virus metagenommics sequencing")
    with report.add_section("Read Statistics", "Read Statistics"):
        df_readcounts = df_clean_read_stats.rename(columns={
            'step': 'Stage',
            'num_seqs': 'Reads',
            'sum_len': 'Total Nucleotides',
            'min_len': 'Shortest Read',
            'avg_len': 'Mean Length',
            'max_len': 'Longest Read'
        })
        DataTable.from_pandas(df_readcounts[['Stage', 'Reads', 'Mean Length']], use_index=False, export=False)
        p(
            """
            Reads left after each filtering step.
            """
        )
    with report.add_section("Consensus Targets", "Consensus Targets"):
        # TODO: use tabs to show another table without description but with more numbers
        # (e.g. N-count, # bases called, )
        df_mapstats = df_mapping_stats.copy()
        df_mapstats['Coverage'] = (
            df_mapstats['CalledNucleobases']
            / (df_mapstats['CalledNucleobases'] + df_mapstats['NCount'])
            * 100
        ).round(0).astype(int)
        df_mapstats.rename(
            inplace=True,
            columns={
                'ReferenceLength': 'Length',
                'NCount': 'Ambiguous positions',
                'NumberOfMappedReads': 'Mapped reads',
                'AverageCoverage': 'Average read coverage',
                'CalledNucleobases': 'Positions called'
            }
        )
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

    # assert len(args.assembly_modes) == len(args.assembly_read_stats)
    # assert len(args.assembly_modes) == len(args.blast_hits)

    os.makedirs(args.outdir, exist_ok=True)
    df_mapping_stats = read_mapping_stats(args.mapping_stats, args.blast_hits)
    df_mapping_stats.to_csv(os.path.join(args.outdir, f'{args.prefix}_consensus_stats.tsv'), sep='\t')
    df_clean_read_stats = pd.read_csv(args.clean_read_stats, sep='\t')

    html_report(
        df_mapping_stats,
        df_clean_read_stats,
        os.path.join(args.outdir, f'{args.prefix}.html')
    )

    # TODO
    # - include assembly report
    # - per filter
    # -- number of reads
    # -- parameters chosen (min_readlen, min_)
    # -- succesful or not
    # -- number of contigs found

    # # TODO: replace this!
    # fnames = [
    #     args.mapping_stats,
    #     args.clean_read_stats,
    #     *args.assembly_read_stats,
    #     *args.blast_hits
    # ]
    # with open(args.report, 'w') as f_report:
    #     f_report.write(" ".join(fnames) + '\n')
    #     for fn in fnames:
    #         with open(fn) as infile:
    #             f_report.write(infile.read() + '\n')


# def main(args):
#     """Run the entry point."""
#     logger = get_named_logger("Report")
#     report = labs.LabsReport(
#         "Influenza Sequencing Report", "wf-flu",
#         args.params, args.versions)
#     with open(args.metadata) as metadata:
#         sample_details = sorted([
#             {
#                 'alias': d['alias'],
#                 'type': d['type'],
#                 'barcode': d['barcode']
#             } for d in json.load(metadata)
#         ], key=lambda d: d["alias"])
#     sample_files = gather_sample_files(sample_details, args)
#     with report.add_section("Typing", "Typing"):
#         p(
#             """
#             This table gives the influenza type and strain for each sample. Samples are
#             first aligned to IRMA to generate a consensus of alignments, and then
#             typed with Abricate using the INSaFLU database. Please see the table in the
#             section below ('Typing details') for full Abricate results. These results
#             are especially useful if typing results are discordant.
#             """
#         )
#         typing_df = typing(sample_details, sample_files)
#         typing_df.columns = typing_df.columns.str.title().str.replace("_", " ")
#         typing_df.rename(columns={'Ha': 'HA', 'Na': 'NA'}, inplace=True)
#         # add a new column 'Archetype'
#         typing_df['Archetype'] = typing_df.apply(lambda row: get_archetype(row), axis=1)
#         typing_df = typing_df[['Alias', 'Barcode', 'Type', 'Archetype']]

#         with table(cls="table"):
#             with thead():
#                 for columns in typing_df.columns:
#                     th(f"{columns}")
#             with tbody():
#                 for index, row in typing_df.iterrows():
#                     with tr():
#                         for i in range(4):
#                             cell = td()
#                             if row[i] == 'undetermined':
#                                 cell.add(h5(span(
#                                     'Undetermined', cls="badge bg-warning"),
#                                     cls="mb-0"))
#                             elif row.index[i] == 'Type' and row[i] == 'Type_A':
#                                 cell.add(
#                                     h5(span('A', cls="badge bg-primary"), cls="mb-0"))
#                             elif row.index[i] == 'Type' and row[i] == 'Type_B':
#                                 cell.add(h5(span('B', cls="badge bg-info"), cls="mb-0"))
#                             elif row.index[i] == 'Archetype' and row.Type == "Type_A":
#                                 cell.add(h5(span(
#                                     row[i], cls="badge bg-primary"), cls="mb-0"))
#                             elif row.index[i] == 'Archetype' and row.Type == "Type_B":
#                                 cell.add(
#                                     h5(span(row[i], cls="badge bg-info"), cls="mb-0"))
#                             elif row.index[i] == 'Barcode' and row[i] is None:
#                                 cell.add(span("None"))
#                             else:
#                                 cell.add(span(row[i]))

#     with report.add_section("Typing details", "Typing details"):
#         tabs = Tabs()
#         with tabs.add_dropdown_menu('Typing details', change_header=True):
#             for sample, files in sample_files.items():
#                 if not os.path.exists(files["type_txt"]):
#                     continue
#                 df = pd.read_csv(files["type_txt"], sep="\t", header=0, index_col=0)
#                 df = df.drop([
#                     'START',
#                     'END',
#                     'STRAND',
#                     'DATABASE',
#                     'ACCESSION',
#                     'PRODUCT'], axis=1)
#                 df.rename(columns={'RESISTANCE': 'DETAILS'}, inplace=True)
#                 df.columns = [x.title() for x in df.columns]
#                 df.columns = df.columns.str.replace('_', ' ')
#                 with tabs.add_dropdown_tab(sample):
#                     DataTable.from_pandas(df, use_index=False, export=True)

#         p(
#             """
#             Select samples from the drop-down in this table to view detailed Abricate
#             results.
#             """
#         )

#     with report.add_section("Coverage", "Coverage"):
#         dfs = []
#         for sample, files in sample_files.items():
#             if not os.path.exists(files["depth"]):
#                 continue
#             df = pd.read_csv(files["depth"], sep="\t", header=None, index_col=0)
#             df.columns = ['position', sample]
#             df.drop(columns='position', inplace=True)
#             df.index.name = 'segment'
#             median = df.groupby('segment').median()
#             dfs.append(median)

#         data = pd.concat(dfs, join='outer', axis=1)
#         data = data.rename_axis("sample", axis=1)

#         plot = ezc.heatmap(data, vmin=0, annot=False)

#         EZChart(plot, 'epi2melabs')

#         p(
#             """
#             The heatmap shows the median coverage per segment for each sample.
#             Each box in the heatmap represents one segment in a sample and is
#             colour-coded using the range of values in the slider (from zero to maximum
#             median coverage across the whole batch). The slider can be manipulated
#             to filter the heatmap by coverage levels, enabling a quick assessment of
#             the coverage for each sample.
#             """
#         )

#     with report.add_section('Nextclade results', 'Nextclade', True):
#         nextclade_data = {}
#         with open(args.nextclade_datasets, "r") as n:
#             csv_reader = csv.DictReader(n, delimiter=",")
#             for record in csv_reader:
#                 nxt_json = f"nextclade/{record['dataset']}.json"
#                 if nxt_json in args.nextclade_files:
#                     output = {
#                         "sample": [],
#                         "strain": [],
#                         "gene": [],
#                         "coverage": [],
#                         "clade": [],
#                         "warnings": []
#                         }
#                     try:
#                         nxt_results = json.load(open(nxt_json))
#                     except json.decoder.JSONDecodeError:
#                         logger.error(
#                             f"Unable to load JSON for {record['dataset']}"
#                         )
#                         sys.exit(1)
#                     for nxt_result in nxt_results["results"]:
#                         output['sample'].append(nxt_result['seqName'])
#                         output['strain'].append(record['strain'])
#                         output['gene'].append(record['gene'])
#                         output['coverage'].append(
#                             f"{nxt_result['coverage']:.2f}"
#                             )
#                         output['clade'].append(nxt_result['clade'])
#                         output['warnings'].append(
#                             ",".join(nxt_result['warnings']))
#                     nextclade_data[f"{record['strain']} - {record['gene']}"] = output
#         if nextclade_data:
#             tabs = Tabs()
#             for tab_name, results in nextclade_data.items():
#                 with tabs.add_tab(tab_name):
#                     df = pd.DataFrame.from_dict(results)
#                     df.columns = df.columns.str.title()
#                     DataTable.from_pandas(
#                         df,
#                         use_index=False,
#                         export=True,
#                         file_name='wf-flu-nextclade')
#         else:
#             p(
#                 """
#                 No sample met the criteria for Nextclade analysis.
#                 Please check the quality and strains of your samples.
#                 Nextclade analysis is available for the following strains:
#                 H1N1pdm, H3N2, Victoria, Yamagata
#                 """
#             )

#     if args.stats:
#         with report.add_section("Read summary", "Read summary"):
#             fastcat.SeqSummary(args.stats)

#     report.write(args.report)
#     logger.info(f"Report written to {args.report}.")

