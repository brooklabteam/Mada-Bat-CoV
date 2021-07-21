# Mada-Bat-CoV

This github repo contains documentation of analyses for the Madagascar bat coronavirus paper. Goal of this repo is to document analysis to define which samples from fecal and urine datasets have detectable coronavirus signal, to build a phylogeny of those CoVs, and to assemble full genomes for as many as possible. Additionally, we will undertake a few simple ecological analyses to highlight the seasonality of bat virus shedding.

1. First, we downloaded all non-host contigs mapping to ANY taxon from IDseq in the 'RR034B1_feces' project, the 'RR034B2_urine_wholeblood_novaseq' project, and the 'RR034B_throat_swab_raw_RNASeq_NovaSeq'project. 

2. Next, we downloaded all the (a) nucleotide and (b) protein full genome **reference sequences** under (i) Betacoronavirus (taxid:694002) and (ii) Alphacoronavirus (taxid: 693996) from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). 

3. Next, we used bowtie2 to map each of the NCBI combined fasta files into a reference dataset.

4. Then, we kicked off a command line blast for each contig subset, blasting these to the reference locally. I ran these on the UC Berkeley computing cluste. Basically, I am re-doing a more focused version of IDseq to see if any other "hits" to CoVs specifically shake out.

