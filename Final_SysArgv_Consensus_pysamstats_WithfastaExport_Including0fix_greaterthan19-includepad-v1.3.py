
# coding: utf-8

# In[1]:

import importlib
import pysam
import pysamstats
import csv
import sys

# run: python scriptname.py filename(without_bam_extention) chromosome fastafile(ref)

# filename of bam file (do not include.bam extention) = sys.argv[1]

# chromosome in file (name after > in bam file) believe this can be removed completely as is only there for situations with multiple chromosomes, we only really have one. = sys.argv[2]

# fastafile = sys.argv[3]
#


# In[3]:

# print(sys.argv[2])
filename = sys.argv[1]
header = ['Position', 'A', 'C', 'G', 'T', 'Insertions', 'Deletions', 'Consensus']
bamfile = pysam.AlignmentFile('{}.bam'.format(filename))
results = []
concensus = []
for record in pysamstats.stat_variation(bamfile,fafile=sys.argv[2],pad=True):
    rec = [record['pos'],record['A'],record['C'],record['G'],record['T'],record['insertions'],record['deletions']]
    if rec[1]+rec[2]+rec[3]+rec[4]:
        percentage = float(max(rec[1:5])) / (rec[1]+rec[2]+rec[3]+rec[4])
    else:
        percentage = 0
        #print(record['pos'])
    ind = rec.index(max(rec[1:5]))
    if ind==1:
        found='A'
    elif ind==2:
        found='C'
    elif ind==3:
        found='G'
    elif ind==4:
        found='T'
    else:
        found='N'
    if max(rec[1:5])>19 and percentage>=0.7:
        rec.append(found)
        concensus.append(found)
    else:
        rec.append('N')
        concensus.append('N')
    results.append(rec)
    seq = ''.join(concensus) # create a string of GTACN

# with open('{}AlignmentFilePileupToRef.csv'.format(filename), 'w', newline='') as csvfile:
# with open('{}AlignmentFilePileupToRef.csv'.format(filename), 'w') as csvfile:
#     writer = csv.writer(csvfile, delimiter=',')
#     writer.writerow(i for i in header)
#     writer.writerows(results)
with open('{}PileupConsensus.fasta'.format(filename), 'w') as fastafile:
    fastafile.write('>{}\n'.format(filename))
    fastafile.write(seq)



