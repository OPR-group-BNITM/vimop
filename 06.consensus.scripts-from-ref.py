import pysam
import pysamstats
import pandas as pd
import csv

# print(snakemake.input[0])
# print(snakemake.input[1])
# print(snakemake.input[2])


# g = open(snakemake.output[3], 'a+')

bamfile = pysam.AlignmentFile(snakemake.input[0])
fastafile = snakemake.input[1]
consensus = []
results = []
# indices_to_loop_through = [1, 2, 3, 4, 6]

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
    if (max(rec[1:5]) > 19 and percentage >= 0.7):
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


seq = ''.join(consensus) # create a string of GTACN
ref_length=len(seq)
percent_ATCG=(nb_ATCG/ref_length)
with open(fastafile, 'r') as f:
    gbtitle = ((f.readline().split(maxsplit=1))[1]).replace("$", "").replace(",", "").replace(";", "").rstrip("\n")   #.rstrip("$")

with open(snakemake.input[2], 'r') as f:
    lines=f.readlines()
    nb_trim_bases=((lines[1].split())[4]).replace(",", "")
    nb_trim_reads=((lines[1].split())[3]).replace(",", "")
    nb_trim_minlen=((lines[1].split())[5]).replace(",", "")
    nb_trim_avglen=((lines[1].split())[6]).replace(",", "")
    nb_trim_maxlen=((lines[1].split())[7]).replace(",", "")


df_clean = pd.read_csv(snakemake.input[3])
data = df_clean['step'].array


cols = pd.MultiIndex.from_product([data, ['Clean']])

df = pd.DataFrame(data, columns=cols)
df['num_seqs'] = df_clean['num_seqs']
# with open(snakemake.input[3], 'r') as f:
#     lines=f.readlines()
#     nb_basesP2=((lines[1].split())[4]).replace(",", "")
#     nb_readsP2=((lines[1].split())[3]).replace(",", "")

print(df)

# with open(snakemake.input[6], 'r') as f:
#     lines=f.readlines()
#     ref_bases=((lines[1].split())[4]).replace(",", "")

# with open(snakemake.input[7], 'r') as f:
#     lines=f.readlines()
#     nb_virus_reads=((lines[4].split())[0]).replace(",", "")

# coverage = pd.read_table(snakemake.input[8], names=['ref','pos','coverage'])
# nb_virus_bases_mapped=coverage['coverage'].sum()

# total_sample_reads= float(nb_readsP1)+float(nb_readsP2)
# fraction_viral_reads=float(nb_virus_reads)/(float(nb_readsP1)+float(nb_readsP2))
# total_sample_bases=float(nb_basesP1)+float(nb_basesP2)
# frac_viral_bases=float(nb_virus_bases_mapped)/(float(nb_basesP1)+float(nb_basesP2))




# with open(snakemake.output[0], 'w') as f:
#     f.write(str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + str(snakemake.params.ref) + "," + str(gbtitle) + "," + str(percent_ATCG) + "," + str(nb_ATCG) + "," + str(ref_length) + "," + str(nb_virus_reads) + "," + str(int(total_sample_reads)) + "," + str(fraction_viral_reads)+ "," + str(int(nb_virus_bases_mapped)) + "," + str(total_sample_bases) + "," + str(frac_viral_bases) + "," + str(seq))

# header = ['Position', 'A', 'C', 'G', 'T', 'Deletions', 'Insertions', 'Consensus']
# with open(snakemake.output[1], 'w') as csvfile:
#     writer = csv.writer(csvfile, delimiter=',')
#     writer.writerow(i for i in header)
#     writer.writerows(results)

# with open(snakemake.output[2], 'w') as f:
#     f.write('>'+ str(snakemake.params.RUNID) + "," + str(snakemake.params.sample) + "," + str(snakemake.params.ref) + '\n' + str(seq))

# f.close

# g.close()
