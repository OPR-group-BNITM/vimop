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
        expand([MAP_PATH + '/{sample}-{target}/02_filter/' + RUNID+'-{sample}-{target}-mapped.fastq', 
            RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-mapped-stats.txt',
            MAP_PATH + '/{sample}-{target}/03_normalize/' + RUNID+'-{sample}-{target}-mapped-normalized.fastq',
            MAP_PATH + '/{sample}-{target}/04_subsample75000/' + RUNID+'-{sample}-{target}-subsampled.fastq'], 
            sample=SAMPLES, target=TARGET)


rule map:
    input:
        DB_DIR+'/{target}.fasta',
        TRIM_PATH+'/'+RUNID+'-{sample}_seqtk_trimfq.fastq'
    output:
        temp(MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.sam')
    benchmark:
        BENCHMARK+'/03_map-{target}/'+RUNID+'-{sample}-{target}-map-benchmark.txt'
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell:
        'minimap2  -a {input[0]} {input[1]} -t {threads} -o {output[0]}'



rule filter:
    input:
        DB_DIR+'/{target}.fasta',
        MAP_PATH + '/{sample}-{target}/01_map/' + RUNID+'-{sample}-{target}-mapped.sam'
    output:
        MAP_PATH + '/{sample}-{target}/02_filter/' + RUNID+'-{sample}-{target}-mapped.fastq',
    conda:
        'envs/general.yaml'
    threads: 4 #workflow.cores
    shell:
        'samtools fastq --threads {threads} -F 4 --reference {input[0]} -1 {output[0]} -2 {output[0]} -n {input[1]} -0 {output[0]} -s {output[0]}'


rule stats_filter:
    input: 
        MAP_PATH + '/{sample}-{target}/02_filter/' + RUNID+'-{sample}-{target}-mapped.fastq'
    output:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-mapped-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input[0]} > {output[0]}'


rule normalize_depth:
    input: 
        MAP_PATH + '/{sample}-{target}/02_filter/' + RUNID+'-{sample}-{target}-mapped.fastq'
    output:
        MAP_PATH + '/{sample}-{target}/03_normalize/' + RUNID+'-{sample}-{target}-mapped-normalized.fastq'
    conda:
        'envs/general.yaml'
    shell:
        'bbnorm.sh in={input[0]} out={output[0]}.fq target=40 mindepth=2'


rule stats_normalize:
    input: 
        MAP_PATH + '/{sample}-{target}/03_normalize/' + RUNID+'-{sample}-{target}-mapped-normalized.fastq'
    output:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-normalized-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input[0]} > {output[0]}'


rule subsample:
    input: 
        MAP_PATH + '/{sample}-{target}/03_normalize/' + RUNID+'-{sample}-{target}-mapped-normalized.fastq'
    output:
        MAP_PATH + '/{sample}-{target}/04_subsample75000/' + RUNID+'-{sample}-{target}-subsampled.fastq'
    conda:
        'envs/general.yaml'
    shell:
        'reformat.sh ow=t samplereadstarget=75000 in={input[0]} out={output[0]}'


rule stats_subsample:
    input: 
        MAP_PATH + '/{sample}-{target}/04_subsample75000/' + RUNID+'-{sample}-{target}-subsampled.fastq'
    output:
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-subsampled-stats.txt'
    conda:
        'envs/general.yaml'
    shell:
        'seqkit stats {input[0]} > {output[0]}'


# rule 

# rule long_read_assembly:
#     input:
#         "long_reads/{sample}_long.fastq.gz"

#     output:
#         "canu-outputs/{sample}.subreads.contigs.fasta"

#     conda:
#         "/home/lamma/env-export/hybrid_assembly.yaml"

#     shell:
#         "canu -p {wildcards.sample} -d canu-outputs genomeSize=8m -pacbio-raw {input}"
# rule short_read_alignment:



# rule spades:
#     input:
#         MAP_PATH + '/{sample}-{target}/02_filter/' + RUNID+'-{sample}-{target}-mapped.fastq',
#     output:
#         touch(MAP_PATH + '/{sample}-{target}/03_assemble/contigs.fasta'),
#         touch(MAP_PATH + '/{sample}-{target}/03_assemble/iflow-spades.log')
#     conda:
#         'envs/general.yaml'
#     params:
#         cov_cutoff=5,
#         kmer_size='21,33,55,77,99,127',
#         outdir=MAP_PATH + '/{sample}-{target}/03_assemble/'

#     benchmark:
#         BENCHMARK+'/03_map-{target}/'+RUNID+'-{sample}-{target}-assemble.txt'
#     threads: 8 #workflow.cores
#     shell:
#         """
#         spades.py \
#         --pe1-1 {input[0]} \
#         --pe1-2 {input[1]} \
#         --cov-cutoff {params.cov_cutoff} \
#         --careful \
#         -t {threads} \
#         -k {params.kmer_size} \
#         -o {params.outdir} > {output[1]} || true
#         """

