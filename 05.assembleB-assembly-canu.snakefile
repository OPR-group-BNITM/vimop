import pandas as pd


RESULTS = config['results']
BENCHMARK = config['benchmark']
TRIM_PATH = config['trim_path']
CLEAN_PATH = config['clean_path']
MAP_PATH = config['map_path']
ASSEMBLE_PATH = config['assemble_path']

RUNID = config['runid']

SAMPLES = (config['currentSample'])

DB_DIR = config['db']


GENOME_SIZE = str(config['genomeSize'])
MIN_READ_LENGTH = str(config['minReadLength'])
MIN_OVERLAP_LENGTH = str(config['minOverlapLength'])
COR_OUT_COVERAGE = str(10000)
ATTEMPT_NUMBER = str(config['attemptNumber'])
SUBSAMPLE_LEVEL = str(config['subsampleLevel'])

rule all:
    input:
        expand([ASSEMBLE_PATH + '/{sample}/02_assemble-canu-'+ATTEMPT_NUMBER+'/' + RUNID+'-{sample}.contigs.fastq'], sample=SAMPLES)



rule canu:
    input:
        ASSEMBLE_PATH + '/{sample}/01_subsample/' + RUNID+'-{sample}-trim-subsampled-'+SUBSAMPLE_LEVEL+'.fastq'
    output:
        touch(ASSEMBLE_PATH + '/{sample}/02_assemble-canu-'+ATTEMPT_NUMBER+'/' + RUNID+'-{sample}.contigs.fastq')
    conda:
        'envs/canu.yaml'
    params:
        outdir = MAP_PATH + '/{sample}/02_assemble-canu-'+ATTEMPT_NUMBER+'/',
        prefix = RUNID + '-{sample}',
        genomeSize = GENOME_SIZE,
        minReadLength = MIN_READ_LENGTH,
        minOverlapLength = MIN_OVERLAP_LENGTH,
        corOutCoverage = COR_OUT_COVERAGE
    benchmark:
        BENCHMARK+'/04_assemble/'+RUNID+'-{sample}-{target}-canu-assemble-'+ATTEMPT_NUMBER+'.txt'
    threads: 8 #workflow.cores
    shell:
        """
        canu \
        -nanopore {input[0]}\
        -fast -d {params.outdir}/ \
        -p {params.prefix} \
        genomeSize={params.genomeSize} \
        minReadLength={params.minReadLength} \
        minOverlapLength={params.minOverlapLength} \
        corOutCoverage={params.corOutCoverage}
        """
