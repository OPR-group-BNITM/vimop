import pandas as pd
from itertools import *
# df=pd.read_csv(config['table'], names = ['sample'])

# SAMPLES = df['sample']
RESULTS = config['results']
BENCHMARK = config['benchmark']
TRIM_PATH = config['trim_path']
CLEAN_PATH = config['clean_path']
MAP_PATH = config['map_path']
ASSEMBLE_PATH = config['assemble_path']
CONSENSUS_PATH = config['consensus_path']

COV_LIMIT = config['covLimit']

RUNID = config['runid']

# TARGET = list((config['target']).split(","))
DB_DIR = config['db']
# ASSEMBLER = list((config['assembler']).split(","))


df=pd.read_csv(config['blastlist'])

STEPS = config['steps'].split(',')
# df1=pd.read_csv(config['all_blasted_list'], names=['sample','ref','def','length','target','assembler','bitscore'])
# df_new=pd.merge(df, df1)

# df_refs=df_new.groupby(['sample','ref']).aggregate({'bitscore': 'sum'}).reset_index()
# df.sort_values(by=['ref'])
# df_refs=df_new[['sample','ref','def','length']]

# df_new = df_new.groupby(['sample','ref','def','length']).aggregate(aggregation_functions).reset_index()
# df_new['def'] = df_new['def'].apply(lambda x: 0 if 'partial' in x else 1)


# df_short = df_new[df_new['length'].apply(lambda x: x in pd.Interval(left=0., right=4000.))]
# df_medium = df_new[df_new['length'].apply(lambda x: x in pd.Interval(left=4001., right=8500.))]
# df_long = df_new[df_new['length'].apply(lambda x: x in pd.Interval(left=8501., right=100000.))]


# df_short = df_short.sort_values(by=['def','bitscore'],ascending=[False, False]).head(n=7)
# df_medium = df_medium.sort_values(by=['def','bitscore'],ascending=[False, False]).head(n=7)
# df_long = df_long.sort_values(by=['def','bitscore'],ascending=[False, False]).head(n=5)

# df_refs = pd.concat([df_short, df_medium, df_long], axis=0, join='outer', ignore_index=False, copy=True)


SAMPLES=df['sample'].tolist()
REF=df['ref'].tolist()

rule all:
  input:
    [CONSENSUS_PATH +'/{sample}/refs/{ref}.fasta'.format(sample=sample_and_ref[0], ref=sample_and_ref[1]) for sample_and_ref in zip(SAMPLES, REF)]
    # [RESULTS +'/{sample}/06_consensus/{RUNID}-{sample}-{ref}-stats.txt'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}-{sample}-{ref}.sam'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}-{sample}-{ref}.bam'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}-sorted-{sample}-{ref}.bam'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}-sorted-{sample}-{ref}.bam.bai'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [RESULTS +'/{sample}/06_consensus/{RUNID}-bam-stats-sorted-{sample}-{ref}.txt'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}-{sample}-{ref}-coverage.txt'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [f'{CONSENSUS_PATH}/{sample_and_ref[0]}/{RUNID}-{sample_and_ref[0]}-{sample_and_ref[1]}-{covlimit}x-consensus.csv' for sample_and_ref, covlimit in product(COV_LIMIT, zip(SAMPLES, REF))],
    # [f'{CONSENSUS_PATH}/{sample_and_ref[0]}/{RUNID}-{sample_and_ref[0]}-{sample_and_ref[1]}-{covlimit}x-alignment-file.csv'for sample_and_ref, covlimit in product(COV_LIMIT, zip(SAMPLES, REF))],
    # [f'{CONSENSUS_PATH}/{sample}/fasta/{RUNID}-{sample_and_ref[0]}-{sample_and_ref[1]}-{covlimit}x-consensus.fasta' for sample_and_ref, covlimit in product(COV_LIMIT, zip(SAMPLES, REF))]
    # [CONSENSUS_PATH +'/{sample}/{RUNID}_{sample}_consensus_all.csv'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}_{sample}_consensus_selected.csv'.format(RUNID=RUNID, sample=sample, ref=ref) for sample, ref in zip(SAMPLES, REF)]



rule get_ref_fasta:
    params:
        virus_db = DB_DIR+'/ALL',
        ref = '{ref}'
    output:
        CONSENSUS_PATH + "/{sample}/refs/{ref}.fasta"
    conda:
        'envs/general.yaml'
    shell: 
        "blastdbcmd -entry {params.ref} -db {params.virus_db} -out {output[0]}"

# rule index_fasta:
#     input:
#         CONSENSUS_PATH + '/{sample}/refs/{ref}.fasta'
#     output:
#         temp(CONSENSUS_PATH + '/{sample}/refs/{ref}.fasta.sa'),
#         temp(CONSENSUS_PATH + '/{sample}/refs/{ref}.fasta.pac'),
#         temp(CONSENSUS_PATH + '/{sample}/refs/{ref}.fasta.bwt'),
#         temp(CONSENSUS_PATH + '/{sample}/refs/{ref}.fasta.ann'),
#         temp(CONSENSUS_PATH + '/{sample}/refs/{ref}.fasta.amb')
#     conda:
#         'envs/general.yaml'
#     shell: 
#         'bwa index {input[0]}'


rule stats_ref_fasta:
    input: 
        CONSENSUS_PATH +'/{sample}/refs/{ref}.fasta'
    output:
        RESULTS+'/{sample}/06_consensus/'+RUNID+'-{sample}-{ref}-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input} > {output}'


rule mapping:
    input: 
        CONSENSUS_PATH +'/{sample}/refs/{ref}.fasta',
        TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq',
    output:
        temp(CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}.sam')
    conda:
        'envs/general.yaml'
    threads: 8 #workflow.cores
    shell:
        'minimap2 -ax map-ont {input[0]} {input[1]} -t {threads} -o {output[0]}'


rule sam_to_bam:
    input: 
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}.sam'
    output:
        temp(CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}.bam')
    conda:
        'envs/general.yaml'
    threads: 1
    shell:
        'samtools view --threads {threads} -b -o {output[0]} {input[0]}'


rule sort_bam:
    input: 
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}.bam'
    output:
        temp(CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam')
    conda:
        'envs/general.yaml'
    threads: 1
    shell:
        'samtools sort --threads {threads} -o {output[0]} {input[0]}'


rule index_sorted_bam:
    input: 
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam'
    output:
        temp(CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam.bai')
    conda:
        'envs/general.yaml'
    shell:
        'samtools index {input[0]}'


rule bam_flagstat:
    input: 
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam'
    output:
        RESULTS+'/{sample}/06_consensus/'+RUNID+'-bam-stats-sorted-{sample}-{ref}.txt'
    conda:
        'envs/general.yaml'
    shell:
        'samtools flagstat {input[0]} > {output[0]}'


rule coverage:
    input: 
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam'
    output:
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-coverage.txt'
    conda:
        'envs/general.yaml'
    shell:
        'bedtools genomecov -ibam {input[0]} -d > {output[0]}'


rule consensus:
    input:
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam',
        CONSENSUS_PATH +'/{sample}/refs/{ref}.fasta',
        RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq-stats.txt',
        RESULTS+'/{sample}/02_clean/' + RUNID+'-{sample}-clean-stats.txt',
        RESULTS+'/{sample}/03_map-'+RUNID+'-{sample}-mapped-assembly-stats.txt',
        RESULTS+'/{sample}/04_assemble/'+RUNID+'-{sample}-canu-assembly-stats.txt',
        RESULTS + '/{sample}/06_consensus/'+RUNID+'-{sample}-{ref}-stats.txt',
        RESULTS + '/{sample}/06_consensus/'+RUNID+'-bam-stats-sorted-{sample}-{ref}.txt',
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-coverage.txt',
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam.bai'
    params:
        ref = '{ref}',
        RUNID = RUNID,
        sample = '{sample}',
        covLimit = COV_LIMIT
    output:
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-{covlimit}x-consensus.csv',
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-{covlimit}x-alignment-file.csv',
        CONSENSUS_PATH +'/{sample}/fasta/'+RUNID+'-{sample}-{ref}-{covlimit}x-consensus.fasta'
        # RESULTS + '/{sample}/{RUNID}-{sample}-{ref}-deletions-positions.csv'
    conda:
        'envs/pysam.yaml'
    script:
        '06.consensus.scripts-from-ref.py'



