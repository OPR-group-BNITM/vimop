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
        expand([RESULTS+'/{sample}/00_initial_stats/'+RUNID+'-{sample}-initial-stats.txt', 
            RESULTS+'/{sample}/00_initial_stats/'+RUNID+'-{sample}-demultiplexed_fastqc.html', 
            TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq-stats.txt',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-seqtk-trimfq_fastqc.html',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.txt',
            RESULTS+'/{sample}/01_trim/'+RUNID+'-{sample}-trimmed-read-length.png'],
            sample=SAMPLES)


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
        DEMULTIPLEX_PATH +'/'+RUNID+'-{sample}-demultiplexed.fastq.gz'
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
        DEMULTIPLEX_PATH +'/'+RUNID+'-{sample}-demultiplexed.fastq.gz'
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
        sample = '{sample}'
    output:
        RESULTS+'/{sample}/01_trim/{RUNID}-{sample}-trimmed-read-length.png'
    script:
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import sys 


        curated_columns=[]
        with open(snakemake.input[0]) as f:
            lines=f.readlines()
            for line in lines:
                if not line.startswith("#"):
                    curated_columns.append(line.split('\t')[0:2])
            


        rl = pd.DataFrame(curated_columns, columns=['Length','Count'])

        # rl = pd.read_csv(snakemake.input[0], sep='\t', names=['Length','Count','Pct_reads','cum_reads', 'cum_pct_reads','bases','pct_bases','cum_bases','cum_pct_bases'])

        # plt.rc('sans-serif':['Arial'], 'size':20})
        fig, ax = plt.subplots(figsize=(20,10))

        fig.set_size_inches(20,10)
        ax.plot(rl, color='xkcd:denim blue',linewidth=3)
        ax.set_facecolor('white')
        # ax.axhline(y=20, linewidth=3, color='black',alpha=0.5, linestyle='dashed', label='20x')

        ax.set_xlabel('Length', size=20,weight='bold')
        ax.set_ylabel('Number of reads', size=20,weight='bold')
        ax.get_xaxis().set_label_coords(0.5,-0.07)
        ax.get_yaxis().set_label_coords(-0.07,0.5)
        ax.tick_params(axis='both', which='major', labelsize=18)

        # ax.title.set_text(str(snakemake.params.RUNID)+"-"+str(snakemake.params.sample)+", reference: "+str(snakemake.params.ref))
        ax.set_title('Run: '+str(snakemake.params.RUNID)+", sample: "+str(snakemake.params.sample)+", trimmed", fontdict={'fontsize': 24, 'fontweight': 'bold'})
        ax.title.set_position((0.5,1.03))
        plt.savefig(snakemake.output[0])


