# Mada-Bat-CoV

This github repo contains documentation of analyses for the Madagascar bat coronavirus paper. Goal of this repo is to document analysis to define which samples from fecal and urine datasets have detectable coronavirus signal, to build a phylogeny of those CoVs, and to assemble full genomes for as many as possible. Additionally, we will undertake a few simple ecological analyses to highlight the seasonality of bat virus shedding.

1. First, we downloaded all non-host contigs mapping to ANY taxon from IDseq in the 'RR034B1_feces' project, the 'RR034B2_urine_wholeblood_novaseq' project, and the 'RR034B_throat_swab_raw_RNASeq_NovaSeq'project. We concatenated each of these into a compiled fasta of non-host reads for each tissue type. Note that these files are too big to include in the GitHub repo.

2. Next, we downloaded all the (a) nucleotide and (b) protein full genome **reference sequences** under (i) Betacoronavirus (taxid:694002) and (ii) Alphacoronavirus (taxid: 693996) from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). We concatenated the AlphaCoV and BetaCoV nucelotide and protein files such that we had only two CoV reference fasta files: one for nt and one for protein. 

3. Next, on the command line, we "cd" into the folder where these are contained and use the following command on each of the fasta files produced in step 2 to build these into two different blast databases: 

```
makeblastdb –in NCBI_all_CoV_nt.fasta –dbtype nucl –parse_seqids
```

Note that you will need to have the command line version of NCBI-Blast installed on your computer to do this. See here for directions.

4. Then, we kicked off a command line blast for each contig subset on the two above databases (so six runs in total). I ran these on the UC Berkeley computing cluste. Basically, I am re-doing a more focused version of IDseq to see if any other "hits" to CoVs specifically shake out.

