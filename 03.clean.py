import os
import sys
import pandas as pd
from shutil import copyfile
from shutil import rmtree

from pathlib import Path
import re

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


for sample in SAMPLES:
    # print(TRIM_PATH+'/'+RUNID+'-'+sample+'_seqtk_trimfq.fastq')
    # print(CLEAN_PATH + '/' + RUNID+'-'+sample+'-pre_'+first_step+'.fastq')
    # os.mkdir(CLEAN_PATH + '/'+sample, mode = 0o777)
    os.makedirs(CLEAN_PATH + '/'+sample, exist_ok=True)
    os.makedirs(ANALYSIS_FOLDER + '/' + RUNID + '_RESULTS/'+sample+'/02_clean/', exist_ok=True)
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
            # os.system('')
for sample in SAMPLES:
    stats = pd.DataFrame(columns=['num_seqs', 'sum_len', 'min_len', 'avg_len','max_len'])
    i = 0
    print(i)
    for step in SPEC2:
        stat_file = open(CLEAN_PATH + '/' + sample+ '/' +step+'/' + RUNID+'-'+sample+'-no-'+step+'-stats.txt','r')

        lines = stat_file.readlines()
        # stats = pd.DataFrame(, data=(re.sub(',','',lines[1])).split()[3:8])
        stats = stats.append(pd.Series((re.sub(',','',lines[1])).split()[3:8], index=['num_seqs', 'sum_len', 'min_len', 'avg_len','max_len']), ignore_index=True
        # print(i)
        # if i == 0:
            # stats['step'] = SPEC[i]
        # else:
            # stats['step'].iloc[-1] = SPEC[i]
        # print(SPEC)
        # print(i)
        # i+=1
    # stats.to_csv(ANALYSIS_FOLDER + '/' + RUNID + '_RESULTS/'+sample+'/02_clean/' + RUNID+'-'+sample+'-clean-stats.txt', index=False)


    # rmtree(CLEAN_PATH + '/' + sample+ '/' +step)


