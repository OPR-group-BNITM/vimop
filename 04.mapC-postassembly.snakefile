import pandas as pd

df=pd.read_csv(config['table'])

SAMPLES = df['sample']
RESULTS = config['results']
BENCHMARK = config['benchmark']
TRIM_PATH = config['trim_path']
CLEAN_PATH = config['clean_path']
MAP_PATH = config['map_path']

RUNID = config['runid']

TARGET = list((config['target']).split(","))
VIRUS_DB = config['virus_db']
ASSEMBLER = config['assembler']



rule stats_assembly:
    input: 
        MAP_PATH + '/{sample}-{target}/05_assemble-{assembler}/contigs.fasta'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-mapped-{assembler}-stats.txt')
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


rule stats_assembly:
    input: 
        MAP_PATH + '/{sample}-{target}/05_assemble-{assembler}/contigs.fasta'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-mapped-{assembler}-stats.txt')
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

rule calc_read_length_trimmed:
    input:
        MAP_PATH + '/{sample}-{target}/05_assemble-{assembler}/contigs.fasta'
    conda:
        'envs/general.yaml'
    output:
        RESULTS +'/{sample}/03_map-{target}/'+RUNID+'-{sample}read-length.txt'
    shell:
        'readlength.sh in={input[0]} in2={input[1]} bin=1 out={output[0]}'


rule trimmed_read_length_png:
    input:
        RESULTS+'/{sample}/ana/read_length/{RUNID}_{sample}_trimmed_read_length.txt'
    conda:
        'envs/general.yaml'
    params:
        RUNID = '{RUNID}',
        sample = '{sample}'
    output:
        RESULTS+'/{sample}/ana/read_length/{RUNID}_{sample}_trimmed_read_length.png'
    script:
        'read_length.py'





rule blast:
    input:
        MAP_PATH + '/{sample}-{target}/03_assemble/contigs.fasta'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-blast-results.lst')
    conda:
        'envs/general.yaml'
    params:
        virus_db=VIRUS_DB+'ALL'
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
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-blast-results.lst'
    output:
        temp(touch(MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_ids'))
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
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-blast-results.lst'
    output:
        temp(touch(MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_bitscore'))
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
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-blast-results.lst'
    output:
        temp(touch(MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_length'))
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
        RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-blast-results.lst'
    output:
        temp(touch(MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_def'))
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
        MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_ids',
        MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_def',
        MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_length',
        MAP_PATH + '/{sample}-{target}/04_blast/'+RUNID+'-{sample}-{target}-hits_bitscore'
    output:
        touch(RESULTS+'/{sample}/03_map-{target}/'+RUNID+'-{sample}-{target}-blasted_list')
    params:
        sample = '{sample}',
        target = '{target}'
    run:
        ref = pd.read_csv(input[0],names=['ref'])
        description = pd.read_csv(input[1],names=['def'])
        length = pd.read_csv(input[2], names=['length'])
        bitscore = pd.read_csv(input[3], names=['bitscore'])
        sam = params.sample
        target = params.target
        merged = pd.concat([ref, description, length, bitscore], axis=1)
        merged['sample']=sam
        merged['target']=target
        merged = merged[['sample','ref','def','length','target','bitscore']]
        aggregation_functions = {'bitscore': 'sum'}
        merged = merged.groupby(['sample','ref','def','length','target']).aggregate(aggregation_functions).reset_index()
        merged.to_csv(output[0], header=False, index=False)



