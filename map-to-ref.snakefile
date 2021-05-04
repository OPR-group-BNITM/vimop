import pandas as pd
from itertools import *


FASTQ = config['fastq']
REFS = config['refs']
OUTDIR = config['outdir']
DB_DIR = config['db']
print(REFS)
print(FASTQ)
# df=pd.read_csv(config['blastlist'], names = ['sample','ref','target'])


# SAMPLES=df['sample'].tolist()
# REF=df['ref'].tolist()

rule all:
  input:
    [OUTDIR +'/{fastq}-{ref}-stats.txt'.format(fastq=fastq, ref=ref) for fastq, ref in product(FASTQ,REFS)]
    # [RESULTS +'/{sample}/05_consensus/{RUNID}-bam-stats-sorted-{sample}-{ref}.txt'.format(RUNID=RUNID, sample=sample_and_ref[0], ref=sample_and_ref[1]) for sample_and_ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH +'/{sample}/{RUNID}-{sample}-{ref}-coverage.txt'.format(RUNID=RUNID, sample=sample_and_ref[0], ref=sample_and_ref[1]) for sample_and_ref in zip(SAMPLES, REF)],
    # [CONSENSUS_PATH+'/{sample}/{RUNID}-{sample}-{ref}-{covlimit}x-consensus.csv'.format(RUNID=RUNID, sample=sample_and_ref[0], ref=sample_and_ref[1], covlimit=covlimit) for sample_and_ref, covlimit in product(zip(SAMPLES, REF), COV_LIMIT)],
    # [CONSENSUS_PATH+'/{sample}/{RUNID}-{sample}-{ref}-{covlimit}x-alignment-file.csv'.format(RUNID=RUNID, sample=sample_and_ref[0], ref=sample_and_ref[1], covlimit=covlimit) for sample_and_ref, covlimit in product(zip(SAMPLES, REF), COV_LIMIT)],
    # [CONSENSUS_PATH+'/{sample}/fasta/{RUNID}-{sample}-{ref}-{covlimit}x-consensus.fasta'.format(RUNID=RUNID, sample=sample_and_ref[0], ref=sample_and_ref[1], covlimit=covlimit) for sample_and_ref, covlimit in product(zip(SAMPLES, REF), COV_LIMIT)]


rule stats_ref_fasta:
    input: 
        OUTDIR +'/{fastq}'
    output:
        OUTDIR +'/{fastq}-{ref}-stats.txt'
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
    threads: 8 
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
        RESULTS+'/{sample}/05_consensus/'+RUNID+'-bam-stats-sorted-{sample}-{ref}.txt'
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
        RESULTS+'/{sample}/03_map-'+RUNID+'-{sample}-canu-mapped-assembly-stats-all-targets.txt',
        RESULTS+'/{sample}/04_assemble/'+RUNID+'-{sample}-canu-assembly-stats.txt',
        RESULTS + '/{sample}/05_consensus/'+RUNID+'-{sample}-{ref}-stats.txt',
        RESULTS + '/{sample}/05_consensus/'+RUNID+'-bam-stats-sorted-{sample}-{ref}.txt',
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-coverage.txt',
        CONSENSUS_PATH +'/{sample}/'+RUNID+'-sorted-{sample}-{ref}.bam.bai'
    params:
        ref = '{ref}',
        RUNID = RUNID,
        sample = '{sample}',
        covlimit = '{covlimit}',
        cleanopts = CLEAN_OPTS

    output:
        touch(CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-{covlimit}x-consensus.csv'),
        touch(CONSENSUS_PATH +'/{sample}/'+RUNID+'-{sample}-{ref}-{covlimit}x-alignment-file.csv'),
        touch(CONSENSUS_PATH +'/{sample}/fasta/'+RUNID+'-{sample}-{ref}-{covlimit}x-consensus.fasta')
    conda:
        'envs/pysam.yaml'
    script:
        '06.consensusB.scripts-from-ref.py'



