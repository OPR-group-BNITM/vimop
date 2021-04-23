import pandas as pd
import glob
import os

df=pd.read_csv(config['table'], names = ['sample'])

SAMPLES = df['sample']
DATA = config['data']
DEMULTIPLEX_PATH = config['demultiplex_path']
TRIM_PATH = config['trim_path']
RESULTS = config['results']
BENCHMARK = config['benchmark']
RUNID = config['runid']

scriptPath = config['script_path']


rule all:
    input:
        expand([TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq-stats.txt',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq_fastqc.html',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.txt',
            # RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.png'
            ],
            sample=SAMPLES)

# rule all:
#     input:
#         expand([RESULTS+'/{sample}/00_initial_stats/'+RUNID+'-{sample}-initial-stats.txt', 
#             RESULTS+'/{sample}/00_initial_stats/'+RUNID+'-{sample}-demultiplexed_fastqc.html', 
#             TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq',
#             RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq-stats.txt',
#             RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq_fastqc.html',
#             RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.txt',
#             RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.png'],
#             sample=SAMPLES)


rule stats_post_demultiplex:
    input:
        DEMULTIPLEX_PATH +'/'+RUNID+'-{sample}-demultiplexed.fastq.gz'
    output:
        RESULTS+'/{sample}/00_initial_stats/'+RUNID+'-{sample}-initial-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input} > {output}'


rule fastqc_stats_post_demultiplex:
    input:
        DEMULTIPLEX_PATH +'/'+RUNID+'-{sample}-demultiplexed.fastq'
    output:
        RESULTS+'/{sample}/00_initial_stats/'+RUNID+'-{sample}-demultiplexed_fastqc.html'
    conda:
        'envs/general.yaml'
    params:
        outdir=RESULTS+'/{sample}/00_initial_stats/'
    threads: 1
    shell:
        'fastqc -q -t {threads} -o {params.outdir} {input}'

    
rule seqtk_trimfq:
    input:
        DEMULTIPLEX_PATH +'/'+RUNID+'-{sample}-demultiplexed.fastq'
    output:
        TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq'
    conda:
        'envs/general.yaml'
    shell:
        'seqtk trimfq -b 30 -e 30 {input} > {output}'


rule stats_post_seqtk_trimfq:
    input:
        TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq'
    output:
        RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input} > {output}'
        

rule fastqc_stats_post_seqtk_trimfq:
    input:
        TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq'
    output:
        RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq_fastqc.html'
    params:
        outdir=RESULTS+'/{sample}/01_trim/'
    conda:
        'envs/general.yaml'
    threads: 1
    shell:
        'fastqc -q -t {threads} -o {params.outdir} {input}'


rule calc_read_length_trimmed:
    input:
        TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq'
    conda:
        'envs/general.yaml'
    output:
        RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.txt'
    shell:
        'readlength.sh in={input[0]} qin=33 bin=1 out={output[0]}'


rule trimmed_read_length_png:
    input:
        RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.txt'
    conda:
        'envs/general.yaml'
    params:
        RUNID = '{RUNID}',
        sample = '{sample}',
        label = 'trimmed'
    output:
        RESULTS+'/{sample}/01_trim/{RUNID}-{sample}-trimmed-read-length.png'
    script:
        scriptPath+'/ana/read-length-individual-image.py'


