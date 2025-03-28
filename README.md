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

## Workflow

In the following, the most important steps of the pipelines are explained with the respective options to set.

### Input

A directory with fastq files is passed to `--fastq`.
The directory can contain subdirectories called barcode01, barcode02, ... 
In this case, the different barcodes are treated as different samples with separate output produced for each.

### Read trimming

At the beginning the pipeline trims the ends of the reads to remove adapter sequences and primers.
`--trim_length` sets the number of bases trimmed from both ends.

### Read classification and optional removal of non-viral reads

Centrifuge is used to classify the reads.
This helps to get an overview of how your metagenomics sample is composed.
Centrifuge also classifies the contigs (see later) to get a rough idea about contigs that do not match anything in the reference data base or only partially. 
Use `--centrifuge_do_classify false` to deactivate all centrifuge classifications and save time.

Additionally, you can use the centrifuge classifications to remove reads (not contigs) that are non-viral to a given degree of confidence.
Activate this with `--centrifuge_do_filter true` and use `--centrifuge_filter_min_score 150` to set the minimum level of confidence (e.g. 150 in our example).
The pipeline will remove all reads that are classified to anything that is not viral with at least this score.
Use this if you have a lot of different contaminations that you want to remove such as bacterial reads from a fecal sample.
But beware that centrifuge is limited in it's accuracy and you may also remove some false positives.

### Host and contaminant removal

Host and contaminant reads are removed by mapping them against reference sequences and extracting those that map.
In our reference data base human DNA and RNA, a mouse genome, a mastomy genome and a set of common lab reagent sequences are included.
The data base config file (contamination.yaml) defines the reference sets and the key values assigned to them.
Use the option `--contamination_filters "reagent,mouse"` for example to remove mouse and reagent reads (the default is humand reads).

### Virus filters

One can also filter for specific species or groups of viruses.
Only reads that map to the given targets are then used in the following assembly step.
An arbitrary number of filters can be used.
The filters themselves are part of the reference data base and the respective names defined in the data base configs.
This command `--targets "MARV,EBOV,FILO"` would activate filters for Marburg virus, Ebola and the Filo-virus family. 

### Assembly

An assembly is run for each of the filtered read sets (explained in the previous section) and for the read set with no target filter.
To disable running the no target filter set `--assemble_notarget false`.

After the initial assembly, an iterative re-assembly procedure is run (unless deactivate with `--reassemble_max_iter 0`).
The purpose is to also find virus segments, that are present in very low concentration.
The reads are then mapped to the contigs and those that map are removed.
The rest is once more assembled.
This procedure is repeated until a maximum number of cycles is reached, no reads are left after filtering or no contigs are produced by the canu assembly.

For read sets that were filtered to a target (see previous section), there is a special procedure if no contigs were assembled.
In this case, the longest X reads (set with `--nocontigs_max_reads_precluster`) are passed to a clustering by cd-hit-est.
Of these clusters, the longest Y reads (set with `--nocontigs_nreads`) are chosen instead of contigs for the following target search.
Canu corrected reads are used if available, else the raw reads.

A number of parameters for the canu assembler can be set. See options.

### Target search

Each contig is used to for a blast search in the virus reference data base.
The highest scoring hit is then used as a reference genome.

### Consensus creation

Reads are mapped against the reference genome.
The mapping parameters can be changed (see Options).
There are three options to generate the consensus.
The default is medaka_variant which uses medaka to call variants and then introduce them to the reference to build the consensus.
The second option is called medaka and uses the medaka consensus functionality.
In this case you don't get a vcf file with variants.
The same is true for the third option called simple using samtools consensus and a correction script for masked areas.

For medaka you can pass a model name for the option `medaka_consensus_model`.
Chosing auto will let medaka chose the model, which is the default.

## Output and report

The output is structured like this:

```
output
├── nf-report.html
├── params.json
└── barcode01
    ├── report_barcode01.html
    ├── classification
    │   ├── classification.html
    │   ├── classification.tsv
    │   ├── classification_kraken.tsv
    │   └── classification_report.tsv
    ├── consensus
    │   ├── AB627954.consensus.fasta
    │   ├── AB627954.depth.txt
    │   ├── AB627954.reads.bam
    │   ├── AB627954.reads.bam.bai
    │   ├── AB627954.reference.fasta
    │   ├── AB627954.variants.vcf.gz
    │   ├── CS272305.consensus.fasta
    │   ├── CS272305.depth.txt
    │   ├── CS272305.reads.bam
    │   ├── CS272305.reads.bam.bai
    │   ├── CS272305.reference.fasta
    │   └── CS272305.variants.vcf.gz
    ├── tables
    │   ├── consensus.tsv
    │   ├── contigs.tsv
    │   └── reads.tsv
    ├── assembly
    │   ├── no-target.contigs.fasta
    │   └── re-assembly.contigs.fasta
    └── selected_consensus
        ├── LCMV_L.fasta
        └── LCMV_S.fasta
```

nf-report.html contains technical information about the run and ressource usage.
For each sample there is a directory with results (here barcode01).
The summary of the results is found in report_samplename.html.
Tables contains the information from the htmml report in .tsv files.
The classification directory contains all files for the centrifuge read classification.
The directory consensus contains consensus genome sequences, alignment files and variants as well as the chosen reference.
The chosen selected genomes for curated virus species are listed in selected_consensus with a separate file for each segment.
Contigs can be found in the fasta files in the assembly directory.

## Options

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

ViMOP relies on a reference data base structure.
It is usually placed in your home directory in a folder called `ViMOP_DB`.
I has the following structure:

```
ViMOP_DB/
├── centrifuge/
├── contaminants/
└── virus/
```

The three subdirectories contain files for centrifuge classification, contaminants/host read removal and the virus reference sequences.
Each directory contains a file with a yaml file with the same name prefix (e.g. centrifuge.yaml, contaminants.yaml, virus.yaml).
The configs hold the relefant information about the data base parts as well as an entry 'version' with a version number and an entry description with a brief 'description'.

The three data base parts are briefly described in the following.

### centrifuge

The centrifuge config looks like this

```yaml
version: 1.0
description: "Refseq reference genomes plus genbank virus sequences"
index_name: all
files:
- all.1.cf
- all.2.cf
- all.3.cf
- all.4.cf
virus_taxid_file: virus_taxids.txt
```

The index name has to be the prefix of the centrifuge index which are the files under files.
These files need to be in the centrifuge DB directory.
The virus_taxid_file contains all virus taxids.
This information is important for the centrifuge based filtering.
In addition to unclassified reads, reads classified to these Tax-IDs will be kept since they are considered to be virus reads.
Version and description are for display in the report.

### contaminants

This directory holds files with sequences of host or reagents.
The respective config file looks like this

```yaml
filters:
  reagent: "reagent-db.fasta.gz"
  human_rna: "GCF_000001405.39_GRCh38.p13_rna.fna.gz"
  human_dna: "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
  mouse: "GCA_000001635.8_GRCm38.p6_genomic.fna.gz"
  mastomys: "GCF_008632895.1_UCSF_Mcou_1_genomic.fna.gz"
version: 1.0
description: "Human (GRCh38), mouse (8_GRCm38), mastomys and contaminant filter set"
```

The keys are used to chose the filters using the command `--contamination_filters`.

### virus

The virus data base contains virus reference sequences.
It consists of a config file, a set of sequence files and a blast data base.
Let's have a look at a small example.
Our directory content could look like this

```
ViMOP_BUCKET/
├── db.yaml
├── ALL.fasta
├── EBOV.fasta
├── LASV.fasta
├── FILO.fasta
└── blast_db/
    ├── ALL.nhr
    ├── ALL.ndb
    ├── ALL.nin
    ├── ALL.njs
    ├── ALL.nog
    ├── ALL.not
    ├── ALL.nos
    ├── ALL.nsq
    ├── ALL.ntf
    └── ALL.nto
```

And the corresponding config could look like this:

```
all:
  blast_db: blast_db
  blast_prefix: ALL
  fasta: ALL.fasta
curated:
  LASV:
    fasta: LASV.fasta
    name: Mammarenavirus lassense
    organisms:
    - Lassa virus GA391
    - Lassa virus Josiah
    - Mammarenavirus lassaense
    segments:
    - S
    - L
  EBOV:
    fasta: EBOV.fasta
    name: Orthoebolavirus
    organisms:
    - Orthoebolavirus bombaliense
    - Orthoebolavirus bundibugyoense
    - Orthoebolavirus restonense
    - Orthoebolavirus sudanense
    - Orthoebolavirus taiense
    - Orthoebolavirus zairense
    segments:
    - Unsegmented
filters:
  FILO: FILO.fasta
version: 1.0
description: "A database that has Lassa and Ebolavirus curated."
params.fasta_sequences: /data/home/nils.petersen/dev/VirusDatasetCuration/workflow/testset/testset.fasta
params.taxa_config: /data/home/nils.petersen/dev/VirusDatasetCuration/workflow/testset/test_groups_refs_and_organisms.yaml
params.filter_max_n_share: 0.01
params.filter_min_relative_length: 0.8
params.cdhit_threshold: 0.98
```

ALL.fasta contains all sequences.
EBOV.fasta contains ebola virus sequences, LASV Lassa virus and FILO filo virus seqeunces.
These are the files that are used as mapping filters.
The file ALL.fasta is also used to create the blast data base.
The curated viruses will be shown in their own sections in the report and have segment wise one fasta file, where the consensus sequence with the highest recovery is chosen.
All other virus targets are added in ALL and 

The headers in the fasta files have to be formatted like in this example

```
>OQ791477.1 |MAG: Lassa mammarenavirus isolate c0212 segment L genomic sequence|Arenaviridae|Mammarenavirus lassaense|forward|L
CAAAATGGGCAACAAGCAAGCCAAGTCAACCAAAGTCGATGAACAACATAGAGCTCAT...
```

Separated with a "|" we have
- genbank ID
- description
- family
- species name
- orientation of the sequence with respect to the original data base entry. We re-oriented sequences so that all sequences of a curated data set have the same orientation. However, this can also simply be set to "Unkown".
- the segment name. Set to "Unknown" for non-curated sequences. For curated sets (e.g. in our example LASV and EBOV) this needs to be assigned. If there is only one segment, use "Unsegmented". The segments also need to be listed in the config file.


### Build your own data base

TODO

## Acknowledgements

This product includes software developed by Oxford Nanopore Technologies Plc.
