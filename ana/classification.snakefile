import pandas as pd

df=pd.read_csv(config['table'])

SAMPLES = df['sample']
# RESULTS = config['results']
# BENCHMARK = config['benchmark']

# TRIM_PATH = config['trim_path']
# CLEAN_PATH = config['clean_path']
DB_DIR = config['db']

RUNID = config['runid']

INPUT = config['input']
OUTDIR = config['outdir']
TAG = config['tag']


rule all:
    input:
        expand([OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-virus.txt',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-report-virus.tsv',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-kraken-style-report-virus.tsv',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-virus.html',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.txt',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-report-all.tsv',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-kraken-style-report-all.tsv',
            OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.html'],
            sample=SAMPLES)


rule classify_virus:
    input:
        INPUT
    output:
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-report-virus.tsv'
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-virus.txt'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/centrifuge.yaml'
    threads: 1 #workflow.cores
    shell:
        'centrifuge -x {params.virus_db}/viral -1 {input[0]} -2 {input[1]} --report-file {output[0]} -S {output[1]}'


rule kraken_style_report_virus:
    input:
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-virus.txt'
    output:
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-kraken-style-report-virus.tsv'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/centrifuge.yaml'
    threads: 1 #workflow.cores
    shell:
        'centrifuge-kreport -x {params.virus_db}/viral {input[0]} > {output[0]}'


rule krona_representation_virus:
    input:
        RESULTS+'/{sample}/ana/classification/'+RUNID+'-{sample}-centrifuge-classification-kraken-style-report-virus.tsv'
    output:
        RESULTS+'/{sample}/ana/classification/'+RUNID+'-{sample}-centrifuge-classification-virus.html'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/general.yaml'
    threads: 1 #workflow.cores
    shell:
        'ktImportTaxonomy -tax $HOME/anaconda3/opt/krona/taxonomy -m 3 -t 5 {input[0]} -o {output[0]}'



rule classify_all:
    input:
        INPUT
    output:
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-report-all.tsv'
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.txt'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/centrifuge.yaml'
    threads: 1 #workflow.cores
    shell:
        'centrifuge -x {params.virus_db}/hpvc -1 {input[0]} -2 {input[1]} --report-file {output[0]} -S {output[1]}'


rule kraken_style_report_all:
    input:
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.txt'
    output:
        OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-kraken-style-report-all.tsv'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/centrifuge.yaml'
    threads: 1 #workflow.cores
    shell:
        'centrifuge-kreport -x {params.virus_db}/hpvc {input[0]} > {output[0]}'


rule krona_representation_all:
    input:
        RESULTS+'/{sample}/ana/classification/'+RUNID+'-{sample}-centrifuge-classification-kraken-style-report-all.tsv'
    output:
        RESULTS+'/{sample}/ana/classification/'+RUNID+'-{sample}-centrifuge-classification-all.html'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/general.yaml'
    threads: 1 #workflow.cores
    shell:
        'ktImportTaxonomy -tax $HOME/anaconda3/opt/krona/taxonomy -m 3 -t 5 {input[0]} -o {output[0]}'



