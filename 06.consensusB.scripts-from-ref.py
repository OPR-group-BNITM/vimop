import pysam
import pysamstats
import pandas as pd
import csv
import os

bamfile = pysam.AlignmentFile(snakemake.input[0])
fastafile = snakemake.input[1]



with open(fastafile, 'r') as f:
    gbtitle = ((f.readline().split(maxsplit=1))[1]).replace("$", "").replace(",", "").replace(";", "").rstrip("\n")   #.rstrip("$")

with open(snakemake.input[2], 'r') as f:
    lines=f.readlines()
    nb_trim_bases=((lines[1].split())[4]).replace(",", "")
    nb_trim_reads=((lines[1].split())[3]).replace(",", "")
    nb_trim_minlen=((lines[1].split())[5]).replace(",", "")
    nb_trim_avglen=((lines[1].split())[6]).replace(",", "")
    nb_trim_maxlen=((lines[1].split())[7]).replace(",", "")


with open(snakemake.input[6], 'r') as f:
    lines=f.readlines()
    ref_bases=((lines[1].split())[4]).replace(",", "")

with open(snakemake.input[7], 'r') as f:
    lines=f.readlines()
    # print(lines[1])
    nb_virus_reads=((lines[4].split())[0]).replace(",", "")

    # nb_virus_bases_mapped=((lines[1].split())[4]).replace(",", "")

coverage = pd.read_table(snakemake.input[8], names=['ref','pos','coverage'])
nb_virus_bases_mapped=coverage['coverage'].sum()

# total_sample_reads= float(nb_readsP1)+float(nb_readsP2)
fraction_viral_reads=float(nb_virus_reads)/(float(nb_trim_reads))

frac_viral_bases=float(nb_virus_bases_mapped)/(float(nb_trim_bases))


if os.path.getsize(snakemake.input[5]) > 0:
    with open(snakemake.input[5], 'r') as f:
        lines=f.readlines()
        nb_assemble_bases=((lines[1].split())[4]).replace(",", "")
        nb_assemble_reads=((lines[1].split())[3]).replace(",", "")
        nb_assemble_minlen=((lines[1].split())[5]).replace(",", "")
        nb_assemble_avglen=((lines[1].split())[6]).replace(",", "")
        nb_assemble_maxlen=((lines[1].split())[7]).replace(",", "")
else:
        nb_assemble_bases=0
        nb_assemble_reads=0
        nb_assemble_minlen=0
        nb_assemble_avglen=0
        nb_assemble_maxlen=0




predf = {
    "RUNID": [snakemake.params.RUNID],
    "Sample": [snakemake.params.sample],
    "Reference": [snakemake.params.ref],
    "NCBI definition": [gbtitle],
    # "Percent ATCG": [percent_ATCG],
    # "Nb base called": [nb_ATCG],
    "Nb bases in reference": [ref_bases],
    "Nb of viral reads": [nb_virus_reads],
    "Sample total bases after trim step": [nb_trim_bases],
    "Sample total reads after trim step": [nb_trim_reads],
    "Fraction viral reads": [fraction_viral_reads],
    "Nb of virus bases": [nb_virus_bases_mapped],
    # "Sample total bases": [nb_trim_bases],
    "Fraction viral bases": [frac_viral_bases],
    "Cleaning options": [snakemake.params.cleanopts],
    'Trim stats, num_seqs': [nb_trim_reads],
    'Trim stats, sum_len': [nb_trim_bases],
    'Trim stats, min_len': [nb_trim_minlen],
    'Trim stats, avg_len': [nb_trim_avglen],
    'Trim stats, max_len': [nb_trim_maxlen],
    'Assembly stats, num_seqs': [nb_assemble_reads],
    'Assembly stats, sum_len': [nb_assemble_bases],
    'Assembly stats, min_len': [nb_assemble_minlen],
    'Assembly stats, avg_len': [nb_assemble_avglen],
    'Assembly stats, max_len': [nb_assemble_maxlen]
}

# predf = dict()
# predf = {
#     "RUNID": [snakemake.params.RUNID],
#     "Sample": [snakemake.params.sample],
#     "Reference": [snakemake.params.ref],
#     "NCBI definition": [gbtitle],
#     # "Percent ATCG": [percent_ATCG],
#     # "Nb base called": [nb_ATCG],
#     "Nb bases in reference": [ref_bases],
#     "Nb of viral reads": [nb_virus_reads],
#     "Sample total bases after trim step": [nb_trim_bases],
#     "Sample total reads after trim step": [nb_trim_reads],
#     "Fraction viral reads": [fraction_viral_reads],
#     "Nb of virus bases": [nb_virus_bases_mapped],
#     # "Sample total bases": [nb_trim_bases],
#     "Fraction viral bases": [frac_viral_bases],
#     "Cleaning options": [snakemake.params.cleanopts],
#     'Trim stats, num_seqs': [nb_trim_reads],
#     'Trim stats, sum_len': [nb_trim_bases],
#     'Trim stats, min_len': [nb_trim_minlen],
#     'Trim stats, avg_len': [nb_trim_avglen],
#     'Trim stats, max_len': [nb_trim_maxlen],
# }



consensus = []
results = []
nb_ATCG = 0
for record in pysamstats.stat_variation(bamfile,fastafile,pad=True):
    rec = [record['pos'],record['A'],record['C'],record['G'],record['T'],record['deletions'],record['insertions']]
    ind = (rec[1:5].index(max(rec[1:5])))+1

    if rec[1]+rec[2]+rec[3]+rec[4] != 0:
        percentage = float(rec[ind]) / (rec[1]+rec[2]+rec[3]+rec[4])
    else:
        percentage = 0
    if ind==1:
        found='A'
    elif ind==2:
        found='C'
    elif ind==3:
        found='G'
    elif ind==4:
        found='T'
    
    # else:
    #     found='N'
    if (max(rec[1:5]) >= int(snakemake.params.covlimit) and percentage >= 0.7):
        rec.append(found)
        consensus.append(found)
        if found != '':
            nb_ATCG +=1
        # else:
        #     g.write(str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + str(snakemake.params.ref) + "," + str(rec[0])+'\n')

    else:
        consensus.append('N')
        rec.append('N')

    results.append(rec)

# seq = ''.join(consensus) # create a string of GTACN
# ref_length=len(seq)
percent_ATCG=(nb_ATCG/int(ref_bases)*100)

key_seq = "Sequence"
key_basecalled = "Nb of bases called"
key_percentconsensuscalled = "% consensus called"
predf[key_seq] = [''.join(consensus)]
predf[key_basecalled] = [nb_ATCG]
predf[key_percentconsensuscalled] = [percent_ATCG]
header = ['Position', 'A', 'C', 'G', 'T', 'Deletions', 'Insertions', 'Consensus']

with open(snakemake.output[1], 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(i for i in header)
    writer.writerows(results)

with open(snakemake.output[2], 'w') as f:
    f.write('>'+ str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + str(snakemake.params.ref) + '\n' + str(''.join(consensus)))





if os.path.exists(snakemake.input[3]) and os.path.getsize(snakemake.input[3]) > 0:

    df_clean = pd.read_csv(snakemake.input[3])
    steps = df_clean['step'].array
    name_cols = []
    for i in range(1,len(steps)+1):
        name_cols.append('Cleaning step '+str(i))
    data = ['num_seqs','sum_len','min_len','avg_len','max_len']
    df_clean['col']= name_cols
    df_clean = df_clean.set_index('col')

    cols = pd.MultiIndex.from_product([['Clean'], name_cols,data], sortorder=None)
    df1 = pd.DataFrame(columns=cols)
    for col in name_cols:
        for dta in data:
    #         print(df_clean.loc[step, dta])
            df1['Clean',col,dta] = pd.Series(df_clean.loc[col,dta])
    df1.columns = [', '.join(col).strip() for col in df1.columns.values]


if os.path.exists(snakemake.input[10]) and os.path.getsize(snakemake.input[10]) > 0:

    df_map_stats = pd.read_csv(snakemake.input[10])

    data = ['num_seqs','sum_len','min_len','avg_len','max_len']

    tgts = df_map_stats['target'].array
    df_map_stats = df_map_stats.set_index('target')

    cols = pd.MultiIndex.from_product([['Mapped reads'],tgts,data], sortorder=None)
    df2 = pd.DataFrame(columns=cols)
    for tgt in tgts:
        for dta in data:
    #         print(df_clean.loc[step, dta])
            df2['Mapped reads',tgt,dta] = pd.Series(df_map_stats.loc[tgt,dta])
    df2.columns = [', '.join(col).strip() for col in df2.columns.values]



if os.path.exists(snakemake.input[4]) and os.path.getsize(snakemake.input[4]) > 0:

    df_map = pd.read_csv(snakemake.input[4])

    data = ['num_seqs','sum_len','min_len','avg_len','max_len']

    tgts = df_map['target'].array
    df_map = df_map.set_index('target')

    cols = pd.MultiIndex.from_product([['Mapped contigs'],tgts,data], sortorder=None)
    df3 = pd.DataFrame(columns=cols)
    for tgt in tgts:
        for dta in data:
    #         print(df_clean.loc[step, dta])
            df3['Mapped contigs',tgt,dta] = pd.Series(df_map.loc[tgt,dta])
    df3.columns = [', '.join(col).strip() for col in df3.columns.values]



# with open(snakemake.input[3], 'r') as f:
#     lines=f.readlines()
#     nb_basesP2=((lines[1].split())[4]).replace(",", "")
#     nb_readsP2=((lines[1].split())[3]).replace(",", "")

# print(df)


df = pd.DataFrame()



df = pd.DataFrame.from_dict(predf)
df1['RUNID'] = snakemake.params.RUNID
df2['RUNID'] = snakemake.params.RUNID
df3['RUNID'] = snakemake.params.RUNID
df1['Sample'] = snakemake.params.sample
df2['Sample'] = snakemake.params.sample
df3['Sample'] = snakemake.params.sample

merged = pd.merge(df,df1,on=['RUNID','Sample'])
merged2 = pd.merge(merged,df2,on=['RUNID','Sample'])
merged3 = pd.merge(merged2,df3,on=['RUNID','Sample'])


merged3.to_csv(snakemake.output[0], index = False)

 # str(int(nb_virus_bases_mapped)) + "," + str(total_sample_bases) + "," + str(frac_viral_bases) + "," + str(seq))


# with open(snakemake.output[0], 'w') as f:
    # f.write(str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + str(snakemake.params.ref) + "," + str(gbtitle) + "," + str(percent_ATCG) + "," + str(nb_ATCG) + "," + str(ref_length) + "," + str(nb_virus_reads) + "," + str(int(total_sample_reads)) + "," + str(fraction_viral_reads)+ "," + str(int(nb_virus_bases_mapped)) + "," + str(total_sample_bases) + "," + str(frac_viral_bases) + "," + str(seq))


# f.close

# g.close()
