import pysam
import pysamstats
import pandas as pd
import csv

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
    nb_virus_reads=((lines[4].split())[0]).replace(",", "")

coverage = pd.read_table(snakemake.input[8], names=['ref','pos','coverage'])
nb_virus_bases_mapped=coverage['coverage'].sum()

# total_sample_reads= float(nb_readsP1)+float(nb_readsP2)
fraction_viral_reads=float(nb_virus_reads)/(float(nb_trim_reads))

frac_viral_bases=float(nb_virus_bases_mapped)/(float(nb_trim_reads))



predf = {
    "RUNID": [snakemake.params.RUNID],
    "Sample": [snakemake.params.sample],
    "ref": [snakemake.params.ref],
    "gbtitle": [gbtitle],
    # "Percent ATCG": [percent_ATCG],
    # "Nb base called": [nb_ATCG],
    "ref length": [ref_bases],
    "nb_virus_reads": [nb_virus_reads],
    "total_sample_reads": [nb_trim_reads],
    "fraction_viral_reads": [fraction_viral_reads],
    "nb_virus_bases_mapped": [nb_virus_bases_mapped],
    "total_sample_bases": [nb_trim_bases],
    "frac_viral_bases": [frac_viral_bases],
    # "seq": [seq]
}




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
percent_ATCG=(nb_ATCG/int(ref_bases))

key_seq = "Seq "+snakemake.params.covlimit+"x"
key_basecalled = "Nb base called"+snakemake.params.covlimit+"x"
key_percentconsensuscalled = "% consensus called"+snakemake.params.covlimit+"x"
predf[key_seq] = [''.join(consensus)]
predf[key_basecalled] = [nb_ATCG]
predf[key_percentconsensuscalled] = [percent_ATCG]
header = ['Position', 'A', 'C', 'G', 'T', 'Deletions', 'Insertions', 'Consensus']
with open(snakemake.output[1], 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(i for i in header)
    writer.writerows(results)

with open(snakemake.output[2], 'w') as f:
    f.write('>'+ str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + ''.join(consensus) + '\n' + str(''.join(consensus)))






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



df_map = pd.read_csv(snakemake.input[4])

data = ['num_seqs','sum_len','min_len','avg_len','max_len']
tgts = df_map['target'].array
df_map = df_map.set_index('target')

cols = pd.MultiIndex.from_product([['Map'],tgts,data], sortorder=None)
df2 = pd.DataFrame(columns=cols)
for tgt in tgts:
    for dta in data:
#         print(df_clean.loc[step, dta])
        df2['Map',tgt,dta] = pd.Series(df_map.loc[tgt,dta])

with open(snakemake.input[5], 'r') as f:
    lines=f.readlines()
    nb_assemble_bases=((lines[1].split())[4]).replace(",", "")
    nb_assemble_reads=((lines[1].split())[3]).replace(",", "")
    nb_assemble_minlen=((lines[1].split())[5]).replace(",", "")
    nb_assemble_avglen=((lines[1].split())[6]).replace(",", "")
    nb_assemble_maxlen=((lines[1].split())[7]).replace(",", "")






# with open(snakemake.input[3], 'r') as f:
#     lines=f.readlines()
#     nb_basesP2=((lines[1].split())[4]).replace(",", "")
#     nb_readsP2=((lines[1].split())[3]).replace(",", "")

# print(df)


df = pd.DataFrame()



df = pd.DataFrame.from_dict(predf)
df1['RUNID'] = snakemake.params.RUNID
df2['RUNID'] = snakemake.params.RUNID
df1['Sample'] = snakemake.params.sample
df2['Sample'] = snakemake.params.sample

merged = pd.merge(df,df1,on=['RUNID','Sample'])
merged2 = pd.merge(merged,df2,on=['RUNID','Sample'])

merged2.to_csv(snakemake.output[0])

 # str(int(nb_virus_bases_mapped)) + "," + str(total_sample_bases) + "," + str(frac_viral_bases) + "," + str(seq))


# with open(snakemake.output[0], 'w') as f:
    # f.write(str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + str(snakemake.params.ref) + "," + str(gbtitle) + "," + str(percent_ATCG) + "," + str(nb_ATCG) + "," + str(ref_length) + "," + str(nb_virus_reads) + "," + str(int(total_sample_reads)) + "," + str(fraction_viral_reads)+ "," + str(int(nb_virus_bases_mapped)) + "," + str(total_sample_bases) + "," + str(frac_viral_bases) + "," + str(seq))


# f.close

# g.close()
