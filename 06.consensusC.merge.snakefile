import sys
import pandas as pd
import pathlib
from itertools import islice  
from pathlib import Path
from shutil import copyfile
# import path
import os

RUNID = config['runid']

df=pd.read_csv(config['table'], names = ['Sample'])


SAMPLES = df['Sample'].tolist()
MAP_PATH = config['map_path']
ASSEMBLE_PATH = config['assemble_path']
CONSENSUS_PATH = config['consensus_path']
RESULTS = config['results']
COV_LIMIT = config['covLimit']
label = config['label']
RUNID_wo_label = config['runidwolabel']
DATE = config['CompletionDay']
COMMONVIRUSES = config['CommonViruses']

# TARGET = list((config['target']).split(","))
# TARGET.append('spades')

predf = dict()

consensusdf = pd.DataFrame()
blastList=pd.read_csv(config['blastlist'], names=['Sample', 'ref','target'])

samples_no_blast = list(set(SAMPLES) - set(blastList['Sample'].tolist()))
# print(samples_no_blast)
# df_no_blast = df[df['Sample'] in samples_no_blast]



for index, row in blastList.iterrows():
    # print(row['Sample'])
    # print(row['ref'])
    for covlimit in COV_LIMIT:
        tmp = pd.read_csv(CONSENSUS_PATH + '/'+ row['Sample'] +'/'+RUNID + '-'+ row['Sample']+'-'+row['ref']+ '-'+ str(covlimit) +'x-consensus.csv')
        tmp['Target'] = row['target']
        tmp['Consensus depth requirement'] = str(covlimit)+'x'
        consensusdf = consensusdf.append(tmp,sort=False,ignore_index=True)
        consensusdf = consensusdf.append(tmp, sort = False)
        # print(tmp.to_string())

# df_all = [ pd.read_csv(config['table'], names = ['sample']) 
for samples_nb in samples_no_blast:
        tmp = pd.read_csv(CONSENSUS_PATH + '/'+ samples_nb +'/'+RUNID + '-'+ samples_nb+'-noblast-consensus.csv')
        tmp['Target'] = ''
        tmp['Consensus depth requirement'] = str(covlimit)+'x'
        consensusdf = consensusdf.append(tmp,sort=False,ignore_index=True)
        consensusdf = consensusdf.append(tmp, sort = False)
consensusdf['Sample'] = consensusdf['Sample'].astype(str)





# merged2 = consensusdf.sort_values(by=['RUNID','Sample','Nb of bases called '+str(max(COV_LIMIT))+'x'],ascending=[True,True,False]) #.reset_index()
merged2 = consensusdf.sort_values(by=['RUNID','Sample','Consensus depth requirement','Nb of bases called',],ascending=[True,True,False,False]) #.reset_index()

merged2['Version'] = config['version']
merged2['Completion date'] = DATE
merged2['Label'] = label
merged2['RUNID+label'] = merged2['RUNID']
        # df['Version'] = config['version']
        # df['Completion date'] = DATE
        # df['Label'] = label
        
merged2['RUNID'] = RUNID_wo_label
merged2['Released?'] = ''
merged2['Analysis comments'] = ''
        # df['Fraction consensus called S'] = df['% consensus called S'].div(100)
        # df['Fraction consensus called L'] = df['% consensus called L'].div(100)
        # df.drop(['index','Fraction viral bases L', 'Fraction viral bases S'], axis=1,inplace = True)
merged2.reset_index(drop=True, inplace=True)
# merged2.to_excel(RESULTS+'/'+RUNID+'-all-consensus-all-targets.xlsx',index=False)











### All samples, all targets, some might be duplicate lines if present for multiple targets
merged2.to_excel(RESULTS+'/'+RUNID+'-all-consensus-all-targets.xlsx',index=False)


df_new = merged2.copy()


df_new['NCBI definition'] = df_new['NCBI definition'].astype(str)
df_new['Partial reference?'] = df_new['NCBI definition'].apply(lambda x: 'partial' if 'partial' in x else 'complete')





with open(COMMONVIRUSES, 'r') as f:
    lines = f.read().splitlines() 

    for line in lines:
        if '#' in line:
            continue
        line = line.replace(': ', ':')
        virus = line.split(':')[0]
        keywords = line.split(':')[1] 

        keywords = keywords.replace(' ,', ',')
        keywords = keywords.replace(', ', ',')
        keywords = keywords.replace(',', '|')
        df = virus + 'df'

        df = df_new[df_new['NCBI definition'].str.contains(keywords)]
        df = df.sort_values(by=['RUNID','Sample','Nb of bases called','Partial reference?'],ascending=[True,True,False,True])


        if 'lassa' in virus and not df.empty:
            ### Select the best S and L segments

            df_short = df[df['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=0., right=4000.))].copy()
            df_long = df[df['Nb bases in reference'].apply(lambda x: x in pd.Interval(left=4001., right=9500.))].copy()

            df_short=df_short.rename(columns={"Reference": "Reference S", "Partial reference?": "Partial reference? S", "NCBI definition": "NCBI definition S","Nb bases in reference": "Nb bases in reference S", "Nb of viral reads": "Nb of viral reads S", "Fraction viral bases": "Fraction viral bases S","Fraction viral reads": "Fraction viral reads S","Nb of virus bases": "Nb of virus bases S","Target": "Target S"})
            df_long=df_long.rename(columns={"Reference": "Reference L", "Partial reference?": "Partial reference? L", "NCBI definition": "NCBI definition L","Nb bases in reference": "Nb bases in reference L", "Nb of viral reads": "Nb of viral reads L", "Fraction viral bases": "Fraction viral bases L","Fraction viral reads": "Fraction viral reads L","Nb of virus bases": "Nb of virus bases L","Target": "Target L"})


            df_short=df_short.rename(columns={"% consensus called":"% consensus called S","Nb of bases called": "Nb of bases called S", "Sequence":"Sequence S"})
            df_long=df_long.rename(columns={"% consensus called":"% consensus called L","Nb of bases called": "Nb of bases called L", "Sequence":"Sequence L"})
            df_short = df_short.sort_values(by=['RUNID','Sample','Consensus depth requirement','Nb of bases called S','Partial reference? S'],ascending=[True,True,False,False,True])
            df_long = df_long.sort_values(by=['RUNID','Sample','Consensus depth requirement','Nb of bases called L','Partial reference? L'],ascending=[True,True,False,False,True])

            df_short = df_short.groupby(['RUNID','Sample','Consensus depth requirement']).first()
            df_short = df_short.reset_index()
            df_long = df_long.groupby(['RUNID','Sample','Consensus depth requirement']).first().reset_index()
            df_long = df_long.reset_index()
            df = pd.merge(df_long,df_short, how='outer').sort_values(by=['RUNID','Sample'],ascending=[True,True])
            df['Fraction consensus called S'] = df['% consensus called S'].div(100)
            df['Fraction consensus called L'] = df['% consensus called L'].div(100)
            df.drop(['index','Fraction viral bases L', 'Fraction viral bases S'], axis=1,inplace = True)
            cols1 = ['RUNID+label','RUNID','Label','Sample','% consensus called S','% consensus called L','Released?','Version','Completion date','Analysis comments',
        'Cleaning options','Sample total reads after trim step','Sample total bases after trim step',
        'Nb of viral reads S','Nb of virus bases S','Fraction viral reads S','Target S','Reference S','NCBI definition S','Partial reference? S',
        'Nb bases in reference S','Nb of bases called S','Fraction consensus called S','Sequence S',
        'Nb of viral reads L','Nb of virus bases L','Fraction viral reads L','Target L','Reference L','NCBI definition L','Partial reference? L',
        'Nb bases in reference L','Nb of bases called L','Fraction consensus called L','Sequence L']

        else:
            df['Fraction consensus called'] = df['% consensus called'].div(100)
            cols1 = ['RUNID+label','RUNID','Label','Sample','% consensus called','Released?','Version','Completion date','Analysis comments',
            'Cleaning options','Sample total reads after trim step','Sample total bases after trim step',
            'Nb of viral reads','Nb of virus bases','Fraction viral reads','Target','Reference','NCBI definition','Partial reference?',
            'Nb bases in reference','Nb of bases called','Fraction consensus called','Sequence']

        # samplesWithConsensus = df['Sample'].tolist()
        # for sample in SAMPLES:
        #     if (sample not in samplesWithConsensus):
        #         predf = {"RUNID": [RUNID],
        #         "Sample": [sample]} 
        #         tmp = pd.DataFrame.from_dict(predf)
        #         df = df.append(tmp, sort = False)



        df.reset_index(drop=True, inplace=True)


        cols2 = df.columns.drop(cols1).tolist()
        cols = cols1 + cols2


        df = df[cols]


        df.sort_values(by=['RUNID','Sample'],ascending=[True,True]).to_excel(RESULTS+'/'+RUNID+'-'+virus+'-selected.xlsx',index=False)

#######################

# for covlimit in COV_LIMIT:
#     df_short_new = df_short[df_short['% consensus called S '+str(covlimit)+'x'] >= 0.7]
#     for sample in SAMPLES:
#         fastaname = RUNID+'-'+sample+'-S-refs_above_0.7-'+str(covlimit)+'x.fasta'
#         with open(RESULTS+'/'+fastaname, 'w') as f:
#         # for index, row in islice(df_short[df_short['Sample']==sample].iterrows(), 1, None):
#             for index, row in df_short_new[df_short_new['Sample']==sample].iterrows():
#                 f.write('>'+RUNID+'-'+sample+'-S-'+str(row['Reference S'])+'\n')
#                 f.write(row['Sequence S '+str(covlimit)+'x']+'\n')
#             f.close()
#     df_short_new.to_excel(RESULTS+'/'+RUNID+'-lassa-S-consensus-comparison-above0.7-'+str(covlimit)+'x.xlsx',index=False)


# for covlimit in COV_LIMIT:
#     df_long_new = df_long[df_long['% consensus called L '+str(covlimit)+'x'] >= 0.7]
#     for sample in SAMPLES:
#         fastaname = RUNID+'-'+sample+'-L-refs_above_0.7-'+str(covlimit)+'x.fasta'
#         with open(RESULTS+'/'+fastaname, 'w') as f:
#         # for index, row in islice(df_long[df_long['Sample']==sample].iterrows(), 1, None):
#             for index, row in df_long_new[df_long_new['Sample']==sample].iterrows():
#                 f.write('>'+RUNID+'-'+sample+'-L-'+str(row['Reference L'])+'\n')
#                 f.write(row['Sequence L '+str(covlimit)+'x']+'\n')
#             f.close()
#     df_long_new.to_excel(RESULTS+'/'+RUNID+'-lassa-L-consensus-comparison-above0.7-'+str(covlimit)+'x.xlsx',index=False)


# df_short = df_short.groupby(['RUNID','Sample']).first()
# df_short = df_short.reset_index()
# df_long = df_long.groupby(['RUNID','Sample']).first().reset_index()
# df_long = df_long.reset_index()


# df_end = pd.concat([df_short, df_long], axis=0, join='outer', ignore_index=False, copy=True).sort_values(by=['RUNID','Sample'],ascending=[True,True]) #.reset_index()

# col = []


# df_end=pd.merge(df_long,df_short, on=mergeOn).sort_values(by=['RUNID','Sample'],ascending=[True,True])
# df.shape()
# df_end[all_cols].to_excel(RESULTS+'/'+RUNID+'-lassa-selected.xlsx',index=False)

# for index, row in df_end.iterrows():
#     for covlimit in COV_LIMIT: 
#         fastaname = str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+str(covlimit)+'x.fasta'
#         with open(RESULTS+'/'+fastaname, 'w') as f:
#             f.write('>'+str(row['RUNID'])+'-'+str(row['Sample'])+'-'+str(row['Segment'])+'-mapto-'+str(row['Reference'])+'\n')
#             f.write(row['Sequence '+str(covlimit)+'x'])
#             f.close()




# quicklook40 = merged2['% consensus called '+str(min(COV_LIMIT))+'x'] >= 0.4
# quicklook50 = merged2['% consensus called '+str(min(COV_LIMIT))+'x'] >= 0.5

# quicklook40.to_excel(RESULTS+'/'+RUNID+'-QUICKLOOK-all-consensus-above-40percent-recovered.xlsx',index=False)
# quicklook50.to_excel(RESULTS+'/'+RUNID+'-QUICKLOOK-all-consensus-above-50percent-recovered.xlsx',index=False)



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







