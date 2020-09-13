import sys
import pandas as pd
import pathlib
from itertools import islice  
from pathlib import Path
from shutil import copyfile
import path
import os

RUNID = config['runid']

df=pd.read_csv(config['table'], names = ['Sample'])

SAMPLES = df['Sample'].tolist()

MAP_PATH = config['map_path']
ASSEMBLE_PATH = config['assemble_path']
CONSENSUS_PATH = config['consensus_path']
RESULTS = config['results']
COV_LIMIT = config['covLimit']

TARGET = list((config['target']).split(","))
TARGET.append('spades')

consensusdf = pd.DataFrame()
blastList=pd.read_csv(config['blastlist'], names=['Sample', 'ref','target'])
for index, row in blastList.iterrows():
    for covlimit in COV_LIMIT:
        tmp = pd.read_csv(CONSENSUS_PATH + '/'+ row['Sample'] +'/'+RUNID + '-'+ row['Sample']+'-'+row['ref']+ '-'+ str(covlimit) +'x-consensus.csv')
        tmp['Target'] = row['target']
        consensusdf = consensusdf.append(tmp)

samplesWithConsensus = consensusdf['Sample'].tolist()
for sample in SAMPLES:
    if (sample not in samplesWithConsensus):
        predf = {"RUNID": RUNID,
        "Sample": sample} 
        tmp = pd.DataFrame.from_dict(predf)
        consensusdf = consensusdf.append(predf)


merged2 = consensusdf.sort_values(by=['RUNID','Sample','Nb of bases called '+str(max(COV_LIMIT))+'x'],ascending=[True,True,False]) #.reset_index()
### All samples, all targets, some might be duplicate lines if present for multiple targets
merged2.to_excel(RESULTS+'/'+RUNID+'-all-consensus-all-targets.xlsx',index=False)



quicklook40 = merged2['% consensus called '+str(min(COV_LIMIT))+'x'] >= 0.4
quicklook50 = merged2['% consensus called '+str(min(COV_LIMIT))+'x'] >= 0.5

quicklook40.to_excel(RESULTS+'/'+RUNID+'-QUICKLOOK-all-consensus-above-40percent-recovered.xlsx',index=False)
quicklook50.to_excel(RESULTS+'/'+RUNID+'-QUICKLOOK-all-consensus-above-50percent-recovered.xlsx',index=False)

### Separate the targets
for i in TARGET:
    tmp = merged2[merged2['Target'] == i]
    tmp.to_excel(RESULTS+'/'+RUNID+'-all-consensus-'+i+'.xlsx',index=False)


### Select the best S and L segments

# aggregation_functions = {'Target': lambda x: ' '.join(x)}
# df_new = merged2.groupby(['RUNID','Sample','Reference','NCBI definition','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']).aggregate(aggregation_functions).reset_index()
df_new = merged2
df_new['Partial reference?'] = df_new['NCBI definition'].apply(lambda x: 'partial' if 'partial' in x else 'complete')

df_new_LASV=df_new[df_new['NCBI definition'].str.contains("Lassa|mammarena")]
df_new_DENV=df_new[df_new['NCBI definition'].str.contains("Dengue")]
df_new_DENV= df_new_DENV.sort_values(by=['RUNID','Sample','Nb of bases called '+str(max(COV_LIMIT))+'x','Partial reference?'],ascending=[True,True,False,True])
df_new_DENV = df_new_DENV.groupby(['RUNID','Sample']).first()
df_new_DENV = df_new_DENV.sort_values(by=['RUNID','Sample'],ascending=[True,True]) #.reset_index()
df_new_DENV.to_excel(RESULTS+'/'+RUNID+'-dengue-selected.xlsx',index=False)


df_short = df_new_LASV[df_new_LASV['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=0., right=4000.))].copy()
df_long = df_new_LASV[df_new_LASV['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=4001., right=9500.))].copy()

# df_short=df_end[df_end['Segment'].values=='S']
# dfL=df_end[df_end['Segment'].values=='L']


df_short=df_short.rename(columns={"Reference": "Reference S", "Partial reference?": "Partial reference? S", "NCBI definition": "NCBI definition S","Nb bases in reference": "Nb bases in reference S", "Nb of viral reads": "Nb of viral reads S", "Fraction viral bases": "Fraction viral bases S","Fraction viral reads": "Fraction viral reads S","Nb of virus bases": "Nb of virus bases S","Target": "Target S"})

df_long=df_long.rename(columns={"Reference": "Reference L", "Partial reference?": "Partial reference? L", "NCBI definition": "NCBI definition L","Nb bases in reference": "Nb bases in reference L", "Nb of viral reads": "Nb of viral reads L", "Fraction viral bases": "Fraction viral bases L","Fraction viral reads": "Fraction viral reads L","Nb of virus bases": "Nb of virus bases L","Target": "Target L"})

colS = []
colL = []
seqS = []
seqL = []
for covlimit in COV_LIMIT:
    df_short=df_short.rename(columns={"% consensus called "+str(covlimit)+"x":"% consensus called S "+str(covlimit)+"x","Nb of bases called "+str(covlimit)+"x": "Nb of bases called S "+str(covlimit)+"x", "Sequence "+str(covlimit)+"x":"Sequence S "+str(covlimit)+"x"})

    df_long=df_long.rename(columns={"% consensus called "+str(covlimit)+"x":"% consensus called L "+str(covlimit)+"x","Nb of bases called "+str(covlimit)+"x": "Nb of bases called L "+str(covlimit)+"x", "Sequence "+str(covlimit)+"x":"Sequence L "+str(covlimit)+"x"})

    colS += ['% consensus called S '+str(covlimit)+'x','Nb of bases called S '+str(covlimit)+'x']
    colL += ['% consensus called L '+str(covlimit)+'x','Nb of bases called L '+str(covlimit)+'x']

    seqS += ['Sequence S '+str(covlimit)+'x']
    seqL += ['Sequence L '+str(covlimit)+'x']

cols_to_order0 = ['RUNID','Sample','Sample total reads','Sample total bases']
cols_to_order1 = ['Nb bases in reference S', 'Nb of viral reads S','Fraction viral reads S', 'Nb of virus bases S','Fraction viral bases S', 'Reference S', 'Partial reference? S', 'NCBI definition S', 'Target S']
cols_to_order2 = ['Nb bases in reference L', 'Nb of viral reads L','Fraction viral reads L', 'Nb of virus bases L','Fraction viral bases L', 'Reference L', 'Partial reference? L', 'NCBI definition L', 'Target L']


cols_to_order_all_but_statsS = cols_to_order0 + colS + cols_to_order1 + seqS 
cols_to_order_all_but_stats = cols_to_order0 + colS + colL + cols_to_order1 + seqS + cols_to_order2 + seqL
mergeOn = cols_to_order0  + (df_short.columns.drop(cols_to_order_all_but_statsS).tolist())
all_cols = cols_to_order_all_but_stats + (df_short.columns.drop(cols_to_order_all_but_statsS).tolist())

df_short = df_short.sort_values(by=['RUNID','Sample','Nb of bases called S '+str(max(COV_LIMIT))+'x','Partial reference? S'],ascending=[True,True,False,True])
df_long = df_long.sort_values(by=['RUNID','Sample','Nb of bases called L '+str(max(COV_LIMIT))+'x','Partial reference? L'],ascending=[True,True,False,True])


for covlimit in COV_LIMIT:
    df_short_new = df_short[df_short['% consensus called S '+str(covlimit)+'x'] >= 0.7]
    for sample in SAMPLES:
        fastaname = RUNID+'-'+sample+'-S-refs_above_0.7-'+str(covlimit)+'x.fasta'
        with open(RESULTS+'/'+fastaname, 'w') as f:
        # for index, row in islice(df_short[df_short['Sample']==sample].iterrows(), 1, None):
            for index, row in df_short_new[df_short_new['Sample']==sample].iterrows():
                f.write('>'+RUNID+'-'+sample+'-S-'+str(row['Reference S'])+'\n')
                f.write(row['Sequence S '+str(covlimit)+'x']+'\n')
            f.close()
    df_short_new.to_excel(RESULTS+'/'+RUNID+'-lassa-S-consensus-comparison-above0.7-'+str(covlimit)+'x.xlsx',index=False)


for covlimit in COV_LIMIT:
    df_long_new = df_long[df_long['% consensus called L '+str(covlimit)+'x'] >= 0.7]
    for sample in SAMPLES:
        fastaname = RUNID+'-'+sample+'-L-refs_above_0.7-'+str(covlimit)+'x.fasta'
        with open(RESULTS+'/'+fastaname, 'w') as f:
        # for index, row in islice(df_long[df_long['Sample']==sample].iterrows(), 1, None):
            for index, row in df_long_new[df_long_new['Sample']==sample].iterrows():
                f.write('>'+RUNID+'-'+sample+'-L-'+str(row['Reference L'])+'\n')
                f.write(row['Sequence L '+str(covlimit)+'x']+'\n')
            f.close()
    df_long_new.to_excel(RESULTS+'/'+RUNID+'-lassa-L-consensus-comparison-above0.7-'+str(covlimit)+'x.xlsx',index=False)


df_short = df_short.groupby(['RUNID','Sample']).first()
df_short = df_short.reset_index()
df_long = df_long.groupby(['RUNID','Sample']).first().reset_index()
df_long = df_long.reset_index()


# df_end = pd.concat([df_short, df_long], axis=0, join='outer', ignore_index=False, copy=True).sort_values(by=['RUNID','Sample'],ascending=[True,True]) #.reset_index()

# col = []


df_end=pd.merge(df_long,df_short, on=mergeOn).sort_values(by=['RUNID','Sample'],ascending=[True,True])
# df.shape()
df_end[all_cols].to_excel(RESULTS+'/'+RUNID+'-lassa-selected.xlsx',index=False)

# for index, row in df_end.iterrows():
#     for covlimit in COV_LIMIT: 
#         fastaname = str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+str(covlimit)+'x.fasta'
#         with open(RESULTS+'/'+fastaname, 'w') as f:
#             f.write('>'+str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+'-mapto-'+str(row['Reference'])+'\n')
#             f.write(row['Sequence '+str(covlimit)+'x'])
#             f.close()







# SAMPLES=df_end['Sample'].tolist()
# REF=df_end['Reference'].tolist()
# SEGMENT=df_end['Segment'].tolist()

# pillar_short = pillar_short.groupby(['RUNID','Sample']).first()
# pillar_long = pillar_long.groupby(['RUNID','Sample']).first()

# df_medium = df_medium.sort_values(by=['def','bitscore'],ascending=[False, False]).head(n=7)


# pillar_end = pd.concat([pillar_short, pillar_long], axis=0, join='outer', ignore_index=False, copy=True).sort_values(by=['RUNID','Sample','Segment'],ascending=[True,True,False]).reset_index()

# df_end = df_end[['RUNID','Sample','Segment','Reference','NCBI definition','Partial reference?','Target','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']]
# pillar_end = pillar_end[['RUNID','Sample','Segment','Reference','NCBI definition','Partial?','Target','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence']]



# df_end['Spades successful?'] = df_end['Sample'].apply(lambda x: "Assembled" if (os.path.exists(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') and (os.path.getsize(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') > 0)) else "failed")
# pillar_end['Spades successful?'] = pillar_end['Sample'].apply(lambda x: "Assembled" if (os.path.exists(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') and (os.path.getsize(ASSEMBLE_PATH+'/'+x+'/contigs.fasta') > 0)) else "failed")


  #  pillar_end[i+'mapping successful?'] = pillar_end['Sample'].apply(lambda x: "Mapped" if (os.path.exists(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') and (os.path.getsize(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') > 0)) else "failed")


# for i in TARGET:
#     if i != 'spades':
#         df_end[i+'mapping successful?'] = df_end['Sample'].apply(lambda x: "Mapped" if (os.path.exists(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') and (os.path.getsize(MAP_PATH+'/'+x+'-'+i+'/03_assemble/contigs.fasta') > 0)) else "failed")
        # if 'meta' not in i:
         



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


# df1=pd.read_csv(RESULTS+'/'+RUNID+'-all-consensus.csv',names=['RUNID','Sample','Reference','NCBI definition','% consensus called','Nb of bases called','Nb bases in reference','Nb of viral reads','Sample total reads','Fraction viral reads','Nb of virus bases','Sample total bases','Fraction viral bases','Sequence'])



# merged=pd.merge(df,df1, how="left", on = ['Sample', 'Reference'])
# merged = merged.drop(['Nb bases in reference_x'], axis=1)
# merged = merged.drop(['NCBI definition_x'], axis=1)
# merged = merged.rename(columns={"NCBI definition_y": "NCBI definition", "Nb bases in reference_y": "Nb bases in reference"})







