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
        expand([MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.fastq', 
            RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-mapped-stats.txt',
            MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled-75000.fastq',
            MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled-50000.fastq',
            RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-subsampled-75000-stats.txt',
            RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-subsampled-50000-stats.txt'], 
            sample=SAMPLES, target=TARGET)



rule map:
    input:
        DB_DIR+'/{target}.fasta',
        TRIM_PATH+'/'+RUNID+'-{sample}-seqtk-trimfq.fastq'
    output:
        temp(MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.sam')
    benchmark:
        BENCHMARK+'/03_map-{target}/'+RUNID+'-{sample}-{target}-map-benchmark.txt'
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell:
        'minimap2 -ax map-ont {input[0]} {input[1]} -t {threads} -o {output[0]}'

rule sambam:
    input:
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.sam',
        DB_DIR+'/{target}.fasta.fai'
    output:
        temp(MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.bam')
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell:
        'samtools view -b -o {output[0]} -t {input[1]} --threads {threads} {input[0]}'


rule filter:
    input:
        DB_DIR+'/{target}.fasta',
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.bam'
    output:
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.fastq',
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell:
        'samtools fastq --threads {threads} -F 4 --reference {input[0]} {input[1]} > {output[0]}'


rule stats_filter:
    input: 
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.fastq'
    output:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-mapped-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input[0]} > {output[0]}'


rule subsample75000:
    input: 
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.fastq'
    output:
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled-75000.fastq'
    conda:
        'envs/general.yaml'
    shell:
        'reformat.sh qin=33 ow=t samplereadstarget=75000 in={input[0]} out={output[0]}'

rule stats_subsample75000:
    input: 
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled-75000.fastq'
    output:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-subsampled-75000-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input[0]} > {output[0]}'




rule subsample50000:
    input: 
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.fastq'
    output:
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled-50000.fastq'
    conda:
        'envs/general.yaml'
    shell:
        'reformat.sh qin=33 ow=t samplereadstarget=50000 in={input[0]} out={output[0]}'


rule stats_subsample50000:
    input: 
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-subsampled-50000.fastq'
    output:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-subsampled-50000-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input[0]} > {output[0]}'


