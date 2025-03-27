# ViMOP

Analysis pipeline for virus metagenomics using nanopore sequencing.

This pipeline is developed by the [Outbreak Preparedness and Response team](https://www.bnitm.de/forschung/forschungsgruppen/pathogen/abt-virologie/laborgruppe-duraffour-pahlmann/team) at the Bernhard Nocht Institute for Tropical medicine.
It is used to analyse nanopore reads from untargeted sequencing of viruses such as Lassa, Ebola or Dengue at various sites.

If you have questions, suggestions, want to contribute or have a specific requirement (e.g. for the license) please feel free to contact us.

## Purpose and limitions

The main purpose of this pipeline is the assembly of virus genomes in metagenomics samples.
We have created a reference data base with our favourite viruses. However, you can also create your own.
For information on data bases read further down.
If you need assistance for setting up a reference data set, please contact us.

The pipeline automatically finds well fitting virus genomes and uses them as references to build reference based consensus genomes.
To build a consensus, [medaka](https://github.com/nanoporetech/medaka) or [samtools consensus](https://www.htslib.org/) are used.
This works well for small and medium size RNA viruses such as Lassa, Dengue, Ebola and many others.
However, for large DNA viruses with repetetive regions such as mpox, this approach will not correctly generate complete genomes.
In any case, we recommend carefully reviewing your output (e.g. the alignment .bam files).

## Prerequisites

This pipeline runs best on a powerful laptop or PC. We recommend at least 30 GB RAM and 16 CPUs.

You can run and install the pipeline from the command line using [nextflow](https://www.nextflow.io/) or from the [EPI2ME desktop](https://nanoporetech.com/software/other/epi2me-desktop-application) application from ONT.
If you are using EPI2ME desktop, nextflow and docker are included in the setup of the software. Else you need to install nextflow and [docker](https://www.docker.com/).

## Installation and operation

You can install this pipeline by cloning this repository, running nextflow or using the EPI2ME desktop.
In any case, you need a reference data base installed.
The pipeline will automatically install all dependencies during the first run using docker, given that you computer is connected to the internet.
After the set up of the reference data base and the dependencies, the pipeline does not require an internet connection. 

### Using the command line

To run the pipeline from the command line, make sure to have nextflow and docker installed.
If you cloned this repository, change into its root directory and run `nextflow main.nf` with the additional options.
Without manually cloning the repository you can simply run `nextflow run OPR-group-BNITM/vimop` plus options.
We will stick to the latter now.
Type `nextflow run OPR-group-BNITM/vimop --set-up-default-db` to install our latest data base release.
It may take a while, since the database has to be downloaded and it's quite huge.
Afterwards, run an analysis with `nextflow run OPR-group-BNITM/vimop --fastq /path/to/fastqfiles --out_dir /path/for/your/output`.
You can get some demo data here: TODO!

### Using EPI2ME desktop

Open the application and enter the github url of this repository under Launch -> Import workflow -> Import from Github.
Once you added the workflow, launch it but check the box Input Options -> Install latest default data base.
This process will download and install the data base into your home directory.
This will probably take a while.
Afterwards you can launch the pipeline to analyse data.
You can also click "run demo analysis".

### Options

The following options can be passed.

| Parameter                             | Default          | Description                                                          |
|---------------------------------------|------------------|----------------------------------------------------------------------|
| `fastq`                               | *(none)*         | Path to input FASTQ files directory                                  |
| `sample_sheet`                        | *(none)*         | A CSV file used to map barcodes to sample aliases                    |
| `out_dir`                             | `output`         | Output directory where sample output subdirectories are created      |
| `output_cleaned_reads`                | `false`          | Output cleaned reads after contamination removal                     |
| `custom_ref_fasta`                    | *(none)*         | Optional custom FASTA file for consensus generation                  |
| `base_db`                             | *                | Path to base directory containing virus and contaminants DB          |
| `virus_db`                            | *                | Path to virus database for BLAST and consensus creation              |
| `virus_db_config`                     | *                | YAML config file for the virus database                              |
| `contaminants_db`                     | *                | Path to contaminants database                                        |
| `contaminants_db_config`              | *                | YAML config for contaminants database                                |
| `classification_db`                   | *                | Path to centrifuge classification DB (.cf files)                     |
| `classification_db_config`            | *                | YAML config for centrifuge classification database                   |
| `trim_length`                         | `30`             | Number of bases trimmed from both ends of each read                  |
| `contamination_filters`               | `reagent,human_dna,human_rna` | Contaminant types to filter                             |
| `centrifuge_do_classify`              | `true`           | Activates the centrifuge read and contig classification              |
| `centrifuge_do_filter`                | `false`          | Enable filtering based on centrifuge classifications                 |
| `centrifuge_filter_min_score`         | `100`            | Minimum classification score to filter non-viral reads               |
| `targets`                             | *(none)*         | Comma-separated target virus list                                    |
| `assemble_notarget`                   | `true`           | Whether to assemble reads not assigned to specific targets           |
| `canu_min_read_length`                | `200`            | Minimum read length for Canu assembly                                |
| `canu_min_overlap_length`             | `200`            | Minimum overlap length for Canu assembly                             |
| `canu_genome_size`                    | `30000`          | Estimated genome size for Canu                                       |
| `canu_read_sampling_bias`             | `5`              | Bias for selecting longer reads during downsampling                  |
| `canu_max_input_coverage`             | `200`            | Max input coverage for Canu assembly                                 |
| `canu_cor_out_coverage`               | `10000`          | Max reads corrected in Canu correction step                          |
| `canu_stop_on_low_coverage`           | `0`              | Stop if coverage drops below this in Canu                            |
| `canu_min_input_coverage`             | `0`              | Min required coverage to run Canu                                    |
| `nocontigs_enable_read_usage`         | `true`           | Use reads even if no contigs assembled                               |
| `nocontigs_max_reads_precluster`      | `10000`          | Max reads for clustering in no-contigs path                          |
| `nocontigs_nreads`                    | `100`            | Number of reads to sample in no-contigs mode                         |
| `nocontigs_cdhit_thresh`              | `0.9`            | CD-HIT threshold for no-contigs                                      |
| `nocontigs_cdhit_wordlen`             | `10`             | CD-HIT word length for no-contigs                                    |
| `consensus_method`                    | `medaka_variant` | Method used for consensus generation                                 |
| `reassemble_cdhit_thresh`             | `0.98`           | CD-HIT threshold for reassembly                                      |
| `reassemble_cdhit_wordlen`            | `10`             | CD-HIT word length for reassembly                                    |
| `reassemble_max_iter`                 | `2`              | Max reassembly iterations                                            |
| `map_to_target_minimap_window_size`   | `5`              | Minimap2 window size for target mapping                              |
| `map_to_target_minimap_kmer_size`     | `11`             | Minimap2 k-mer size for seed generation                              |
| `medaka_consensus_model`              | `auto`           | Medaka model name (overrides auto-selection)                         |
| `consensus_min_depth`                 | `20`             | Minimum coverage depth for consensus base                            |
| `consensus_min_share`                 | `0.7`            | Minimum allele share for consensus base                              |
| `min_ram_gb`                          | `30`             | Minimum RAM for workflow                                             |
| `min_cpus`                            | `16`             | Minimum number of CPUs for workflow                                  |
| `min_disk_space_work_gb`              | `100`            | Min disk space for working directory                                 |
| `min_disk_space_out_gb`               | `10`             | Min disk space for output directory                                  |

\* All these files have defaults corresponding to the default data base. If one wants to replace them, one can pass a different path here.

## Database

TODO

## Acknowledgements

This product includes software developed by Oxford Nanopore Technologies Plc.
