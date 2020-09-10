import sys
import pandas as pd
import pathlib
from itertools import islice  
from pathlib import Path
from shutil import copyfile
import path
import os

RUNID = config['runid']

df=pd.read_csv(config['table'])

SAMPLES = df['sample'].tolist()

MAP_PATH = config['map_path']
ASSEMBLE_PATH = config['assemble_path']
CONSENSUS_PATH = config['consensus_path']
RESULTS = config['results']
COV_LIMIT = config['covLimit']

TARGET = list((config['target']).split(","))
# for target in TARGET:
#     for j in SAMPLES:
#         if os.path.isfile(MAP_PATH+'/'+j+'-'+target+'/03_assemble/contigs.fasta'):
#             copyfile(MAP_PATH+'/'+j+'-'+target+'/03_assemble/contigs.fasta', RESULTS+'/'+RUNID+'-'+j+'-mapped-'+target+'-contigs.fasta')
#         if os.path.isfile(ASSEMBLE_PATH+'/'+j+'/contigs.fasta'):
#             copyfile(ASSEMBLE_PATH+'/'+j+'/contigs.fasta', RESULTS+'/'+RUNID+'-'+j+'-assembled-contigs.fasta')

TARGET.append('spades')

# grep -hs ^ $ANALYSIS_FOLDER/${RUNID}-05-consensus/*/${RUNID}_*_consensus.csv \
# > $ANALYSIS_FOLDER/${RUNID}_RESULTS/${RUNID}-all-consensus.csv


consensusdf = pd.DataFrame()
blastList=pd.read_csv(config['all_blasted_list'])
for index, rows in blastList.iterrows():
    for covlimit in COV_LIMIT:
        tmp = pd.read_csv(CONSENSUS_PATH + '/'+ row['sample'] +'/'+RUNID + '-'+ row['sample']+'-'+row['ref']+ '-'+ covlimit +'x-consensus.csv')
        consensusdf.append(tmp)
for sample in SAMPLES:
    if (sample not in consensusdf['Sample'].tolist()):
        predf = {"RUNID": RUNID,
        "Sample": sample} 
        tmp = pd.DataFrame.from_dict(predf)
        consensusdf.append(predf)

# df1=pd.read_csv(RESULTS+'/'+RUNID+'-all-consensus.csv',names=['RUNID','Sample','Reference','NCBI definition','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence'])



# merged=pd.merge(df,df1, how="left", on = ['Sample', 'Reference'])
# merged = merged.drop(['Nb bases in reference_x'], axis=1)
# merged = merged.drop(['NCBI definition_x'], axis=1)
# merged = merged.rename(columns={"NCBI definition_y": "NCBI definition", "Nb bases in reference_y": "Nb bases in reference"})



merged2 = consensusdf.sort_values(by=['RUNID','Sample','Nb of bases called'],ascending=[True,True,False]).reset_index()
### All samples, all targets, some might be duplicate lines if present for multiple targets
merged2.to_excel(RESULTS+'/'+RUNID+'-all-consensus-all-targets.xlsx',index=False)

quicklook40 = merged2['% consensus called'] >= 0.4
quicklook50 = merged2['% consensus called'] >= 0.5

quicklook40.to_excel(RESULTS+'/'+RUNID+'-QUICKLOOK-all-consensus-above-40percent-recovered.xlsx',index=False)
quicklook50.to_excel(RESULTS+'/'+RUNID+'-QUICKLOOK-all-consensus-above-50percent-recovered.xlsx',index=False)

### Separate the targets
for i in TARGET:
    tmp = merged2[merged2['Target'] == i]
    tmp.to_excel(RESULTS+'/'+RUNID+'-all-consensus-'+i+'.xlsx',index=False)
    tmp = merged2[merged2['Target'] == i+'-meta']
    tmp.to_excel(RESULTS+'/'+RUNID+'-all-consensus-'+i+'-meta.xlsx',index=False)




### Select the best S and L segments

aggregation_functions = {'Target': lambda x: ' '.join(x)}
df_new = merged2.groupby(['RUNID','Sample','Reference','NCBI definition','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']).aggregate(aggregation_functions).reset_index()
df_new['Partial reference?'] = df_new['NCBI definition'].apply(lambda x: 'partial' if 'partial' in x else 'complete')

df_new_LASV=df_new[df_new['NCBI definition'].str.contains("Lassa|mammarena")]
df_new_DENV=df_new[df_new['NCBI definition'].str.contains("Dengue")]
df_new_DENV= df_new_DENV.sort_values(by=['RUNID','Sample','Nb of bases called','Partial reference?'],ascending=[True,True,False,True])
df_new_DENV = df_new_DENV.groupby(['RUNID','Sample']).first()
df_new_DENV = df_new_DENV.sort_values(by=['RUNID','Sample'],ascending=[True,True]).reset_index()
df_new_DENV.to_excel(RESULTS+'/'+RUNID+'-dengue-selected.xlsx',index=False)


# df_new_pillar = merged #['RUNID','Sample','Reference','NCBI definition','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']
# df_new_pillar=df_new_pillar[df_new_pillar['NCBI definition'].str.contains("Lassa|mammarena")]
# df_new_pillar = df_new_pillar[~df_new_pillar['Target'].str.contains('meta')]
# df_new_pillar = df_new_pillar.groupby(['RUNID','Sample','Reference','NCBI definition','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']).aggregate(aggregation_functions).reset_index()
# df_new_pillar['Partial?'] = df_new_pillar['NCBI definition'].apply(lambda x: "partial" if "partial" in x else "complete")


df_short = df_new_LASV[df_new_LASV['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=0., right=4000.))].copy()
df_long = df_new_LASV[df_new_LASV['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=4001., right=9500.))].copy()
df_short['Segment']='S'
df_long['Segment']='L'


# for index, row in df_short.iterrows():
#     fastaname = str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+'-mapto-'+str(row['Reference'])+'.fasta'
#     with open(RESULTS+'/'+row['Sample']+'/05_consensus/'+fastaname, 'w') as f:
#         f.write('>'+fastaname+'\n')
#         f.write(row['Sequence'])
#         f.close()


# for index, row in df_long.iterrows():
#     fastaname = str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+'-mapto-'+str(row['Reference'])+'.fasta'
#     with open(RESULTS+'/'+row['Sample']+'/05_consensus/'+fastaname, 'w') as f:
#         f.write('>'+fastaname+'\n')
#         f.write(row['Sequence'])
#         f.close()


# pillar_short = df_new_pillar[df_new_pillar['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=0., right=4000.))].copy()
# pillar_long = df_new_pillar[df_new_pillar['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=4001., right=9500.))].copy()
# pillar_short['Segment']='S'
# pillar_long['Segment']='L'


# pillar_short = pillar_short[~df_short['Target'].str.contains('meta')]


df_short = df_short.sort_values(by=['RUNID','Sample','Nb of bases called','Partial reference?'],ascending=[True,True,False,True])
df_long = df_long.sort_values(by=['RUNID','Sample','Nb of bases called','Partial reference?'],ascending=[True,True,False,True])

# pillar_short = pillar_short.sort_values(by=['RUNID','Sample','Partial?','Nb of bases called'],ascending=[True,True,True,False])
# pillar_long = pillar_long.sort_values(by=['RUNID','Sample','Partial?','Nb of bases called'],ascending=[True,True,True,False])


df_short_new = df_short[df_short['% consensus called'] >= 0.7]
df_long_new = df_long[df_long['% consensus called'] >= 0.7]


for sample in SAMPLES:
    fastaname = RUNID+'-'+sample+'-S-refs_above_0.7.fasta'
    with open(RESULTS+'/'+fastaname, 'w') as f:
        # for index, row in islice(df_short[df_short['Sample']==sample].iterrows(), 1, None):
        for index, row in df_short_new[df_short_new['Sample']==sample].iterrows():
            f.write('>'+RUNID+'-'+sample+'-S-'+str(row['Reference'])+'\n')
            f.write(row['Sequence']+'\n')
        f.close()


for sample in SAMPLES:
    fastaname = RUNID+'-'+sample+'-L-refs_above_0.7.fasta'
    with open(RESULTS+'/'+fastaname, 'w') as f:
        # for index, row in islice(df_long[df_long['Sample']==sample].iterrows(), 1, None):
        for index, row in df_long_new[df_long_new['Sample']==sample].iterrows():
            f.write('>'+RUNID+'-'+sample+'-L-'+str(row['Reference'])+'\n')
            f.write(row['Sequence']+'\n')
        f.close()


df_short_new.to_excel(RESULTS+'/'+RUNID+'-lassa-S-consensus-comparison-above0.7.xlsx',index=False)
df_long_new.to_excel(RESULTS+'/'+RUNID+'-lassa-L-consensus-comparison-above0.7.xlsx',index=False)


df_short = df_short.groupby(['RUNID','Sample']).first()
df_long = df_long.groupby(['RUNID','Sample']).first()

# pillar_short = pillar_short.groupby(['RUNID','Sample']).first()
# pillar_long = pillar_long.groupby(['RUNID','Sample']).first()

# df_medium = df_medium.sort_values(by=['def','bitscore'],ascending=[False, False]).head(n=7)
df_end = pd.concat([df_short, df_long], axis=0, join='outer', ignore_index=False, copy=True).sort_values(by=['RUNID','Sample','Segment'],ascending=[True,True,False]).reset_index()
# pillar_end = pd.concat([pillar_short, pillar_long], axis=0, join='outer', ignore_index=False, copy=True).sort_values(by=['RUNID','Sample','Segment'],ascending=[True,True,False]).reset_index()

df_end = df_end[['RUNID','Sample','Segment','Reference','NCBI definition','Partial reference?','Target','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']]
# pillar_end = pillar_end[['RUNID','Sample','Segment','Reference','NCBI definition','Partial?','Target','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']]



df_end['Spades successful?'] = df_end['Sample'].apply(lambda x: "Assembled" if (os.path.exists(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') and (os.path.getsize(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') > 0)) else "failed")
# pillar_end['Spades successful?'] = pillar_end['Sample'].apply(lambda x: "Assembled" if (os.path.exists(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') and (os.path.getsize(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') > 0)) else "failed")




for i in TARGET:
    if i != 'spades':
        df_end[i+'mapping successful?'] = df_end['Sample'].apply(lambda x: "Mapped" if (os.path.exists(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') and (os.path.getsize(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') > 0)) else "failed")
        # if 'meta' not in i:
           #  pillar_end[i+'mapping successful?'] = pillar_end['Sample'].apply(lambda x: "Mapped" if (os.path.exists(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') and (os.path.getsize(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') > 0)) else "failed")


df_end.to_excel(RESULTS+'/'+RUNID+'-lassa-selected.xlsx',index=False)

for index, row in df_end.iterrows():
    fastaname = str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+'.fasta'
    with open(RESULTS+'/'+fastaname, 'w') as f:
        f.write('>'+str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+'-mapto-'+str(row['Reference'])+'\n')
        f.write(row['Sequence'])
        f.close()


SAMPLES=df_end['Sample'].tolist()
REF=df_end['Reference'].tolist()
SEGMENT=df_end['Segment'].tolist()














