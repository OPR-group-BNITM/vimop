import pandas as pd

df=pd.read_csv(config['table'], names = ['sample'])


SAMPLES = df['sample']
RESULTS = config['results']
BENCHMARK = config['benchmark']
TRIM_PATH = config['trim_path']
CLEAN_PATH = config['clean_path']
MAP_PATH = config['map_path']

RUNID = config['runid']

TARGET = list((config['target']).split(","))
DB_DIR = config['db']




rule all:
    input:
        expand([MAP_PATH + '/{sample}-{target}/02_assemble-canu/contigs.fasta'], 
            sample=SAMPLES, target=TARGET)



rule canu:
    input:
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled.fastq'
    output:
        touch(MAP_PATH + '/{sample}-{target}/02_assemble-canu/contigs.fasta'),
    conda:
        'envs/canu.yaml'
    params:
        outdir = MAP_PATH + '/{sample}-{target}/02_assemble-canu/',
        prefix = RUNID + '-{sample}-{target}',
        cov_cutoff = 5,
        genomeSize = 10000,
        minReadLength = 300,
        minOverlapLength = 200,
        corOutCoverage = 10000
    benchmark:
        BENCHMARK+'/03_map-{target}/'+RUNID+'-{sample}-{target}-canu-assemble.txt'
    threads: 8 #workflow.cores
    shell:
        """
        canu \
        -nanopore-raw {input[0]}\
        -fast -d {params.outdir}/ \
        -p {params.prefix} \
        -t {threads} \
        genomeSize={params.genomeSize} \
        minReadLength={params.minReadLength} \
        minOverlapLength={params.minOverlapLength} \
        corOutCoverage={params.corOutCoverage}
        """
