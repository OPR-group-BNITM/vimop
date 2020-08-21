import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys 


curated_columns=[]
with open(snakemake.input[0]) as f:
    lines=f.readlines()
    for line in lines:
        if not line.startswith("#"):
            curated_columns.append(line.split('\t')[0:2])
    


rl = pd.DataFrame(curated_columns, columns=['Length','Count'])

# rl = pd.read_csv(snakemake.input[0], sep='\t', names=['Length','Count','Pct_reads','cum_reads', 'cum_pct_reads','bases','pct_bases','cum_bases','cum_pct_bases'])

# plt.rc('sans-serif':['Arial'], 'size':20})
fig, ax = plt.subplots(figsize=(20,10))

fig.set_size_inches(20,10)
ax.plot(rl, color='xkcd:denim blue',linewidth=3)
ax.set_facecolor('white')
# ax.axhline(y=20, linewidth=3, color='black',alpha=0.5, linestyle='dashed', label='20x')

ax.set_xlabel('Length', size=20,weight='bold')
ax.set_ylabel('Number of reads', size=20,weight='bold')
ax.get_xaxis().set_label_coords(0.5,-0.07)
ax.get_yaxis().set_label_coords(-0.07,0.5)
ax.tick_params(axis='both', which='major', labelsize=18)

# ax.title.set_text(str(snakemake.params.RUNID)+"-"+str(snakemake.params.sample)+", reference: "+str(snakemake.params.ref))
ax.set_title('Run: '+str(snakemake.params.RUNID)+", sample: "+str(snakemake.params.sample)+", trimmed", fontdict={'fontsize': 24, 'fontweight': 'bold'})
ax.title.set_position((0.5,1.03))
plt.savefig(snakemake.output[0])