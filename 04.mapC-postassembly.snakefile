import pandas as pd
import os

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
ASSEMBLER = list((config['assembler']).split(","))

# SAMPLES2 = df['sample'].tolist()
# for sample in SAMPLES:
#     for target in TARGET:
#         for assembler in ASSEMBLER:
#     #     for assembler 
#             file=RESULTS +'/'+sample+'/03_map-'+target+'/'+RUNID+'-'+sample+'-'+target+'-'+assembler+'-'+'contigs.fasta'
#             if (os.path.exists(file) and os.path.getsize(file) > 0):
#                 SAMPLES2.remove(sample)

rule all:
    input:
        expand([RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-mapped-assembly-stats.txt',
            RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-read-length.txt',
            # RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-read-length.png',
            RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-results.lst',
            RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blasted-list.csv'], 
            sample=SAMPLES, target=TARGET, assembler= ASSEMBLER)




rule rename:
    input: 
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs.fasta'
    output:
        touch(RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs-renamed.fasta')
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail;
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            seqkit rename  {input[0]} > {output[0]}
        fi
        exit 0
        """





rule sort:
    input: 
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs-renamed.fasta'
    output:
        touch(RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs-sorted.fasta')
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail;
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            seqkit sort --by-length --reverse {input[0]} > {output[0]}
        fi
        exit 0
        """



rule stats_assembly:
    input: 
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs-sorted.fasta'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-mapped-assembly-stats.txt')
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail;
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            seqkit stats {input[0]} > {output[0]}
        fi
        exit 0
        """


rule calc_read_length:
    input:
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs-sorted.fasta'
    conda:
        'envs/general.yaml'
    output:
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-read-length.txt'
    shell:
        'readlength.sh qin=33 in={input[0]} bin=1 out={output[0]}'


# rule read_length_png:
#     input:
#         RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-read-length.txt'
#     conda:
#         'envs/general.yaml'
#     params:
#         RUNID = RUNID,
#         sample = '{sample}',
#         label = 'mapped-{target}-{assembler}'
#     output:
#         RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-read-length.png'
#     script:
#         'ana/read-length-individual-image.py'





rule blast:
    input:
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-contigs-sorted.fasta'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-results.lst')
    conda:
        'envs/general.yaml'
    params:
        virus_db=DB_DIR+'/ALL'
    threads: 4 #workflow.cores
    shell:
        """
        set +o pipefail; 
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            blastn -num_threads {threads} -db {params.virus_db} -query {input[0]} -out {output[0]} -outfmt 5
        fi
        exit 0
        """

 
rule get_hit_id:
    input:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-results.lst'
    output:
        temp(touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-ids'))
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail; 
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            grep -A 8 "<Hit_num>1</Hit_num>" {input[0]} | grep Hit_accession | sed "s/  <Hit_accession>//g" | sed "s/<\/Hit_accession>//g"  > {output[0]}
        fi
        exit 0
        """

rule get_hit_bitscore:
    input:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-results.lst'
    output:
        temp(touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-bitscore'))
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail;
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            grep -A 8 "<Hit_num>1</Hit_num>" {input[0]} | grep bit-score | sed "s/      <Hsp_bit-score>//g" | sed "s/<\/Hsp_bit-score>//g" > {output[0]}
        fi
        exit 0
        """


rule get_hit_length:
    input:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-results.lst'
    output:
        temp(touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-length'))
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail;
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            grep -A 8 "<Hit_num>1</Hit_num>" {input[0]} | grep Hit_len | sed "s/  <Hit_len>//g" | sed "s/<\/Hit_len>//g" > {output[0]}
        fi
        exit 0
        """

rule get_hit_def:
    input:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-results.lst'
    output:
        temp(touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-def'))
    conda:
        'envs/general.yaml'
    shell:
        """
        set +o pipefail;
        if [ -f "{input[0]}" ] && [ -s "{input[0]}" ]; then
            grep -A 8 "<Hit_num>1</Hit_num>" {input[0]} | grep Hit_def | sed "s/  <Hit_def>//g" | sed "s/<\/Hit_def>//g" | sed "s/,//g" | sed "s/;//g" > {output[0]}
        fi
        exit 0
        """


rule merge_id_bitscore:
    input:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-ids',
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-def',
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-length',
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blast-hits-bitscore'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-{assembler}-blasted-list.csv')
    params:
        sample = '{sample}',
        target = '{target}',
        assembler = '{assembler}'
    run:
        ref = pd.read_csv(input[0],names=['ref'])
        description = pd.read_csv(input[1],names=['def'])
        length = pd.read_csv(input[2], names=['length'])
        bitscore = pd.read_csv(input[3], names=['bitscore'])
        sam = params.sample
        target = params.target
        assembler = params.assembler
        merged = pd.concat([ref, description, length, bitscore], axis=1)
        merged['sample']=sam
        merged['target']=target
        merged['assembler']=assembler
        merged = merged[['sample','ref','def','length','target','assembler','bitscore']]
        aggregation_functions = {'bitscore': 'sum'}
        merged = merged.groupby(['sample','ref','def','length','target','assembler']).aggregate(aggregation_functions).reset_index()
        merged.to_csv(output[0], header=False, index=False)

