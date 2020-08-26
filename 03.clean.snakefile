import pandas as pd
from shutil import copyfile

df=pd.read_csv(config['table'], names = ['sample'])

SAMPLES = df['sample']
RESULTS = config['results']
BENCHMARK = config['benchmark']
TRIM_PATH = config['trim_path']
CLEAN_PATH = config['clean_path']

RUNID = config['runid']

DB_DIR = config['db']
step = config['step']


# print(DB_DIR)
MASTER_DB = {'reagent':DB_DIR+'/reagent-db.fasta.gz','human_rna':DB_DIR+'/GCF_000001405.39_GRCh38.p13_rna.fna.gz', 
      'human_dna':DB_DIR+'/GCF_000001405.39_GRCh38.p13_genomic.fna.gz','mouse':DB_DIR+'/GCA_000001635.8_GRCm38.p6_genomic.fna.gz',
     'mastomys':DB_DIR+'/GCF_008632895.1_UCSF_Mcou_1_genomic.fna.gz'}

MAP_READS_CMD = 'minimap2 -x map-ont \
-a \
{params.actual_db} \
{input[0]} \
-t {threads} \
-o {output[0]}'

# FILTER_CMD = 'samtools view --threads {threads} \
# -b -q 50 -U  {output[0]} -o {output[1]} {input[0]}'

FILTER_CMD = 'samtools \
fastq -f 4 \
--reference {params.actual_db} \
-1 {output[0]} -2 {output[0]} -0 {output[0]} -s {output[0]} -n \
{input[0]}'



# FASTQ_CMD = 'samtools \
# fastq -n --threads {threads} {input[0]} \
# -1 {output[0]} -2 {output[1]}'
rule all:
    input:
        expand([CLEAN_PATH + '/{sample}/'+step+'/'+RUNID+'-{sample}-no-'+step+'.fastq'], sample=SAMPLES)

rule remove_db:
    input:
        CLEAN_PATH + '/{sample}/' + RUNID+'-{sample}-pre_'+step+'.fastq'
    params:
        actual_db = MASTER_DB.get(step.split('-')[-1])
    output:
        temp(CLEAN_PATH + '/{sample}/'+step+'/' + RUNID+'-{sample}-'+step+'.aln.sam')
    benchmark:
        BENCHMARK+'/02_clean/'+RUNID+'-{sample}-'+step+'-map-benchmark.txt'
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell: 
        (MAP_READS_CMD)


rule filter_unmapped:
    input:
        CLEAN_PATH + '/{sample}/'+step+'/' + RUNID+'-{sample}-'+step+'.aln.sam'
    params:
        actual_db = MASTER_DB.get(step.split('-')[-1])
    output:
        CLEAN_PATH + '/{sample}/'+step+'/' + RUNID+'-{sample}-no-'+step+'.fastq'
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell:
        FILTER_CMD
