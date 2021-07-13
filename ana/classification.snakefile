import pandas as pd

df=pd.read_csv(config['table'], names = ['sample'])


SAMPLES = df['sample']
# RESULTS = config['results']
# BENCHMARK = config['benchmark']

# TRIM_PATH = config['trim_path']
# CLEAN_PATH = config['clean_path']
DB_DIR = config['db']

RUNID = config['runid']

# INPUT = config['input']
# OUTDIR = config['outdir']
# TAG = config['tag']
TRIM_PATH = config['trim_path']
RESULTS = config['results']

rule all:
    input:
        expand([RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-virus.txt',
            RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-report-virus.tsv',
            RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-kraken-style-report-virus.tsv',
            RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-virus.html',
            # OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.txt',
            # OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-report-all.tsv',
            # OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-kraken-style-report-all.tsv',
            # sOUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.html'
            ],
            sample=SAMPLES)


rule classify_virus:
    input:
        TRIM_PATH +'/'+RUNID+'-{sample}-seqtk-trimfq.fastq'
    output:
        RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-report-virus.tsv',
        RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-virus.txt'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/centrifuge.yaml'
    threads: 1 #workflow.cores
    shell:
        'centrifuge -x {params.virus_db}/viral -U {input[0]} --report-file {output[0]} -S {output[1]}'


rule kraken_style_report_virus:
    input:
        RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-virus.txt'
    output:
        RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-kraken-style-report-virus.tsv'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/centrifuge.yaml'
    threads: 1 #workflow.cores
    shell:
        'centrifuge-kreport -x {params.virus_db}/viral {input[0]} > {output[0]}'


rule krona_representation_virus:
    input:
        RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-kraken-style-report-virus.tsv'
    output:
        RESULTS+'/classification/{sample}/'+RUNID+'-{sample}-centrifuge-classification-virus.html'
    params:
        virus_db=DB_DIR
    conda:
        '../envs/general.yaml'
    threads: 1 #workflow.cores
    shell:
        'ktImportTaxonomy -tax $HOME/opt/krona/taxonomy -m 3 -t 5 {input[0]} -o {output[0]}'



# rule classify_all:
#     input:
#         INPUT
#     output:
#         OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-report-all.tsv',
#         OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.txt'
#     params:
#         virus_db=DB_DIR
#     conda:
#         '../envs/centrifuge.yaml'
#     threads: 1 #workflow.cores
#     shell:
#         'centrifuge -x {params.virus_db}/hpvc -U {input[0]} --report-file {output[0]} -S {output[1]}'


# rule kraken_style_report_all:
#     input:
#         OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.txt'
#     output:
#         OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-kraken-style-report-all.tsv'
#     params:
#         virus_db=DB_DIR
#     conda:
#         '../envs/centrifuge.yaml'
#     threads: 1 #workflow.cores
#     shell:
#         'centrifuge-kreport -x {params.virus_db}/hpvc {input[0]} > {output[0]}'


# rule krona_representation_all:
#     input:
#         OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-kraken-style-report-all.tsv'
#     output:
#         OUTDIR+'/'+RUNID+'-{sample}-'+TAG+'-centrifuge-classification-all.html'
#     params:
#         virus_db=DB_DIR
#     conda:
#         '../envs/general.yaml'
#     threads: 1 #workflow.cores
#     shell:
#         'ktImportTaxonomy -tax $HOME/anaconda3/opt/krona/taxonomy -m 3 -t 5 {input[0]} -o {output[0]}'



