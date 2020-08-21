import os
import sys
import pandas as pd
from shutil import copyfile
from pathlib import Path


ANALYSIS_FOLDER = sys.argv[1]

df=pd.read_csv(sys.argv[1]+'/table.csv', names = ['sample'])

SAMPLES = df['sample']

SPEC = sys.argv[2].split(",")

RUNID = sys.argv[3]
TRIM_PATH = sys.argv[4]
CLEAN_PATH = sys.argv[5]
scriptPath = sys.argv[6]
threads = sys.argv[7]
HOME = str(Path.home())

SPEC2 = []
SPEC2.append(SPEC[0])
for i in range(1,len(SPEC)):
    SPEC2.append(SPEC2[i-1]+'-'+SPEC[i])

first_step = SPEC2[0]
last_step = SPEC2[-1]

# os.mkdir(CLEAN_PATH, mode = 0o777)

for sample in SAMPLES:
    # print(TRIM_PATH+'/'+RUNID+'-'+sample+'_seqtk_trimfq.fastq')
    # print(CLEAN_PATH + '/' + RUNID+'-'+sample+'-pre_'+first_step+'.fastq')
    # os.mkdir(CLEAN_PATH + '/'+sample, mode = 0o777)
    os.makedirs(CLEAN_PATH + '/'+sample, exist_ok=True)
    copyfile(TRIM_PATH+'/'+RUNID+'-'+sample+'-seqtk-trimfq.fastq', CLEAN_PATH + '/'+sample+'/' + RUNID+'-'+sample+'-pre_'+first_step+'.fastq')


for step in SPEC2:
    for sample in SAMPLES:
        os.makedirs(CLEAN_PATH + '/'+sample+'/'+step, exist_ok=True)
    os.system('snakemake \
    --snakefile '+scriptPath+'/03.clean.snakefile \
    --config step='+step+' \
    --configfile '+ANALYSIS_FOLDER+'/config.yaml --use-conda --conda-prefix '+HOME+'/opt/iflow \
    --cores '+ threads)
    for sample in SAMPLES:
        if (step != last_step):
            next_step=SPEC2[(SPEC2.index(step)) + 1]
            copyfile(CLEAN_PATH + '/' + sample+ '/' +step+'/' + RUNID+'-'+sample+'-no-'+step+'.fastq', CLEAN_PATH + '/' + sample+ '/' + RUNID+'-'+sample+'-pre_'+next_step+'.fastq')
        else:
            copyfile(CLEAN_PATH + '/' + sample+ '/' +step+'/' + RUNID+'-'+sample+'-no-'+step+'.fastq', CLEAN_PATH + '/' + sample+ '/' + RUNID+'-'+sample+'-cleaned.fastq')


