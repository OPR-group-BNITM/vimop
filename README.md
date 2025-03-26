# ViMOP

Analysis pipeline for virus metagenomics using nanopore sequencing (e.g. with the sispa-protocol).

## 

TODO:
- initial stuff
- you need to install the software
- set up the data base

## Requirements and Prerequesites/Dependencies

This pipeline runs best on a powerful laptop or PC. We recommend at least 30 GB RAM and 16 CPUs.

You can run and install the pipeline from the command line using nextflow or from the EPI2ME desktop application from ONT (https://nanoporetech.com/software/other/epi2me-desktop-application).
If you are using EPI2ME desktop, nextflow and docker are included in the setup of the software. Else you need to install nextflow and docker.

## Installation

You can install this pipeline by cloning this repository, running nextflow or using the EPI2ME desktop.
If nextflow and docker are installed, you can run this pipeline by typing `nextflow OPR-group-BNITM/`.

TODO: with set up DB

## Run nanoflow 

Todo

## Options

The following options can be passed.

| Parameter                             | Default          | Description                                                                                   |
|---------------------------------------|------------------|-----------------------------------------------------------------------------------------------|
| `fastq`                               | *(none)*         | Path to input FASTQ files directory                                                           |
| `sample_sheet`                        | *(none)*         | A CSV file used to map barcodes to sample aliases                                             |
| `out_dir`                             | `output`         | Output directory where sample output subdirectories are created                               |
| `output_cleaned_reads`                | `false`          | Output cleaned reads after contamination removal                                              |
| `custom_ref_fasta`                    | *(none)*         | Optional custom FASTA file for consensus generation                                           |
| `centrifuge_classification_library`   | `all`            | Centrifuge classification library to use                                                      |
| `base_db`                             | *                | Path to base directory containing virus and contaminants DB                                   |
| `virus_db`                            | *                | Path to virus database for BLAST and consensus creation                                       |
| `virus_db_config`                     | *                | YAML config file for the virus database                                                       |
| `contaminants_db`                     | *                | Path to contaminants database                                                                 |
| `contaminants_db_config`              | *                | YAML config for contaminants database                                                         |
| `classification_db`                   | *                | Path to centrifuge classification DB (.cf files)                                              |
| `trim_length`                         | `30`             | Number of bases trimmed from both ends of each read                                           |
| `contamination_filters`               | `reagent,human_dna,human_rna` | Contaminant types to filter                                                      |
| `classification_virus_taxids`         | *(none)*         | TaxID file to remove non-viral reads                                                          |
| `centrifuge_filter_do_it`             | `false`          | Enable filtering based on centrifuge classifications                                          |
| `centrifuge_filter_min_score`         | `100`            | Minimum classification score to filter non-viral reads                                        |
| `targets`                             | *(none)*         | Comma-separated target virus list                                                             |
| `assemble_notarget`                   | `true`           | Whether to assemble reads not assigned to specific targets                                    |
| `canu_min_read_length`                | `200`            | Minimum read length for Canu assembly                                                         |
| `canu_min_overlap_length`             | `200`            | Minimum overlap length for Canu assembly                                                      |
| `canu_genome_size`                    | `30000`          | Estimated genome size for Canu                                                                |
| `canu_read_sampling_bias`             | `5`              | Bias for selecting longer reads during downsampling                                           |
| `canu_max_input_coverage`             | `200`            | Max input coverage for Canu assembly                                                          |
| `canu_cor_out_coverage`               | `10000`          | Max reads corrected in Canu correction step                                                   |
| `canu_stop_on_low_coverage`           | `0`              | Stop if coverage drops below this in Canu                                                     |
| `canu_min_input_coverage`             | `0`              | Min required coverage to run Canu                                                             |
| `nocontigs_enable_read_usage`         | `true`           | Use reads even if no contigs assembled                                                        |
| `nocontigs_max_reads_precluster`      | `10000`          | Max reads for clustering in no-contigs path                                                   |
| `nocontigs_nreads`                    | `100`            | Number of reads to sample in no-contigs mode                                                  |
| `nocontigs_cdhit_thresh`              | `0.9`            | CD-HIT threshold for no-contigs                                                               |
| `nocontigs_cdhit_wordlen`             | `10`             | CD-HIT word length for no-contigs                                                             |
| `consensus_method`                    | `medaka_variant` | Method used for consensus generation                                                          |
| `reassemble_cdhit_thresh`             | `0.98`           | CD-HIT threshold for reassembly                                                               |
| `reassemble_cdhit_wordlen`            | `10`             | CD-HIT word length for reassembly                                                             |
| `reassemble_max_iter`                 | `2`              | Max reassembly iterations                                                                     |
| `map_to_target_minimap_window_size`   | `5`              | Minimap2 window size for target mapping                                                       |
| `map_to_target_minimap_kmer_size`     | `11`             | Minimap2 k-mer size for seed generation                                                       |
| `medaka_consensus_model`              | `auto`           | Medaka model name (overrides auto-selection)                                                  |
| `consensus_min_depth`                 | `20`             | Minimum coverage depth for consensus base                                                     |
| `consensus_min_share`                 | `0.7`            | Minimum allele share for consensus base                                                       |
| `min_ram_gb`                          | `30`             | Minimum RAM for workflow                                                                      |
| `min_cpus`                            | `16`             | Minimum number of CPUs for workflow                                                           |
| `min_disk_space_work_gb`              | `100`            | Min disk space for working directory                                                          |
| `min_disk_space_out_gb`               | `10`             | Min disk space for output directory                                                           |


## Database

TODO

## Acknowledgements

This product includes software developed by Oxford Nanopore Technologies Plc.