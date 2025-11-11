# Output interpretation
In this tutorial we will go through our demo results and learn how to interpret the output of ViMOP. We will also explain which metrics to look at to evaluate the analysis quality and possibly identify false positives in the results. For this we will go through the HTML report provided by ViMOP section by section.  

After completing the installation tutorial, you will end up with a [report file](https://opr.bnitm.de/example_data/report_example_data.html) and the output files. Feel free to explore our example report through the previously provided link.  
You will find this HTML-report when selecting the **Reports tab** in the EPI2ME interface, or in the output directory that you would have defined via the command line.  

# Table of Contents

1. [Output interpretation](#output-interpretation)
2. [Summary of how to evaluate your results](#summary-of-how-to-evaluate-your-results)
3. [Step by step through the report](#step-by-step-through-the-report)
   - [Read Statistics](#read-statistics)
   - [Consensus](#consensus)
   - [Contigs](#contigs)
4. [How to review unexpected results](#how-to-review-unexpected-results)
   - [Open your output folder](#open-your-output-folder)

# Summary of how to evaluate your results

As an overview: when evaluating your analysis results you should look out for the following cases in your report. We will go through these cases and where to detect them in the upcoming sections.
- **Partial matches** are often indicated by:
  - short **BLAST hit length**
  - low **Reference Alignment coverage**
  - low **Contig alignment coverage**
- **Wrong Reference Selection** could have happened when you see
  - a low **Coverage** of the consensus sequence
  - low **Sequence identity** between contig and BLAST hit
  - **mismatch between centrifuge and BLAST** contig classification
  - The Contig **Length** is much shorter or longer than the BLAST **Hit Length**

In these cases you should manually review the output files of your analysis run. A first step for such a reviewing process will be described in the last section of this tutorial.

# Step by step through the report

## Read Statistics
Once you open the report, the first thing you will see are the read statistics.  
This includes:
- **Read length distribtions** before and after trimming and filtering 
  - check, if the mean read length of your sample was enough for a reliable consensus
  - ViMOP can handle mean read lengths as low as 100 bp, but more ideal would be any value >500 bp.
- Amount of **filtered reads** based on the centrifuge classification and the mapping onto the host genome and the reagent database
- **Read classification plot**
  - here we can see the metagenomic overview of our sample
  - If canu failed to build a contig and thus the workflow did not yield a result, you can check if there was any virus found here and restart the analysis with a custom fasta file as a reference. 


Our demo run only contains 5000 simulated reads, a real data set would have much more reads, but the length and quality distribution would be similar since the simulator was trained on a real LASV sequencing data set.  
In this section you can   
In the Read classification tab we can see the metagenomic overview of our sample. If canu failed to build a contig and thus the workflow did not yield a result, you can check if there was any virus found here and restart the analysis with a custom fasta file as a reference. 

![image.png](03_example_output_interpretation_files/image.png)  

![image-2.png](03_example_output_interpretation_files/image-2.png)  

![image-3.png](03_example_output_interpretation_files/image-3.png)  

![image-4.png](03_example_output_interpretation_files/image-4.png)


## Consensus
This section reports the main output of the pipeline. If ViMOP found a reference in the target detection step, the best BLAST hit per virus species would be used for consensus sequence creation.  
We have a curated subset of viruses relevant to our projects. If these viruses are found in our samples, these viruses show up in the report in a separate tab.  
In our demo report for example we have both LASV segments and the MS2 phage as a control in the sample. All other viruses would show up in the Non-curated and the All tab.  
To get an idea of how successful the consensus creation was, have a look at the following sections in this table:
- **Overview** tab->**Coverage**
  - indicates how much of the reference could be reconstructed from the sequencing run supported by a minimum read depth of 20X.
  - If the Coverage is very low, the resulting consensus sequence may not be sufficient for phylogenetic analysis and should be used with caution.
- **Details** tab->**Positions called**, **Ambigious positions**, **mapped reads**, **read coverage**
  - these metrics can be indicators for a reliable consensus sequence. The more positions are called with a high average read coverage, the more reliable the consensus sequence
  - if these metrics are low, this is usually an indicator that there were not enough reads in the sample for a reliable consensus sequence. 
  - in the demo report you can see that 997 and 999 reads were used to assemble the LASV-S and LASV-L sequences. This makes sense, because we simulated the data to contain 1000 reads of each LASV segments. This means that our mapping used almost every read that was contained in the sample to build the consensus as good as possible.

![image-5.png](03_example_output_interpretation_files/image-5.png)  
![image-8.png](03_example_output_interpretation_files/image-8.png)

## Contigs
In the **Contigs section** you will see the assembly results provided by canu. You can also see the BLAST hit of each contig that was used to create the consensus sequence. Additionally, there is the centrifuge classification of each contig.  
To review your contigs look at the following metrics:
- **Reference Alignment coverage**
  - if this value is below 1 it means that the contig is only a partial match and does not fully cover the found BLAST hit.
  - a very low value means that canu probably only managed to assemble a fragment of the found genome
- **Contig alignment coverage**
  - a value below 1 it means that the reference is only a partial match to the contig
  - a very low value should raise eyebrows: how much longer is the contig than the BLAST hit?
  - if the contig is significantly longer than the BLAST hit it could be from a bacterial genome that also contains viral genetic material through endogenization and was thus misclassified. This often happens with bacteriophages and their respective prokaryotes. If you see something like this, the first check would be to look up the **centrifuge classification** of this contig
- **Sequence identity**
  - a low sequence identity means that the contig and found BLAST hit do not actually match
  - this should prompt the user to re-evaluate the found reference
- **BLAST Hit** and **centrifuge classification** 
  - check the organism column -- does it contain the same organism as the Classification column?
  - if not we probably have a false positive result of a contig that was misclassified as a virus

![image-6.png](03_example_output_interpretation_files/image-6.png)  

In the **Assembly statistics** tab, you can see how many reads were given to, corrected and then used for assembly by canu. Canu subsamples a small fraction of reads for assembly, since using all of them would take too long. The number of 1534 that we see here is what you would also see for similar runs with more reads.

![image-7.png](03_example_output_interpretation_files/image-7.png)

# How to review unexpected results

## Open your output folder

If you saw any unexpected outputs as described previously, or if ViMOP failed to create a consensus sequence in general, you should:
1. Open the assembly subfolder in your output folder.  
![image-11.png](03_example_output_interpretation_files/image-11.png)
2. Open the FASTA files in that subfolder
3. Select the contig(s) that caused your suspicion and upload its sequence to [NCBI BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch) and BLAST it against the most recent BLAST database version.  
    - If this yields the same (weird) results as ViMOP or no hits at all, then the sequencing run probably did not contain enough viral material for virus detection or Canu failed to create useful contigs.
    - If this search, yields a different virus that fits more into your expectations with a high identity, then the sample contained a virus we proabably have not in our database yet.
4. You could then restart the analysis by providing the found fasta File as a custom reference. 
    - In EPI2ME that would be under Input Options->Custom reference FASTA. 
    - Via the command line you can provide a custom reference FASTA (with one or multiple references) using the --custom**ref**fasta parameter.




