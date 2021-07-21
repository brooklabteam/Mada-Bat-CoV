# Mada-Bat-CoV

This github repo contains documentation of analyses for the Madagascar bat coronavirus paper. Goal of this repo is to document analysis to define which samples from fecal and urine datasets have detectable coronavirus signal, to build a phylogeny of those CoVs, and to assemble full genomes for as many as possible. Additionally, we will undertake a few simple ecological analyses to highlight the seasonality of bat virus shedding.

---

1. First, we downloaded all non-host contigs derived from fecal, urine, or saliva and mapping to ANY taxon from IDseq in the 'RR034B1_feces' project, the 'RR034B2_urine_wholeblood_novaseq' project, and the 'RR034B_throat_swab_raw_RNASeq_NovaSeq'project. (Note that we did not take contigs from HeLa controls or water, and in the case of the 'RR034B2_urine_wholeblood_novaseq' project, we only looked at urine samples). We concatenated each of these into a compiled fasta of non-host reads for each tissue type. Note that these files are too big to include in the GitHub repo.

---

2. Now, because many of the contigs will overlap in sequences and slow our mapping down, we can deduplicate the compiled .fasta for each of the tissues using the program [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/). You'll need to [install this](https://github.com/weizhongli/cdhit/wiki/2.-Installation) on your home computer first. If using MacOS, I recommend using bioconda to do it--see [here](https://anaconda.org/bioconda/cd-hit). Once installed, you can deduplicate each of the contig files via the following command line script in the same folder as the downloaded contigs (here shown just for the fecal contigs file -- you will need to run it for all three tissue types). This command :

```
cd-hit-est -i rr034b1_feces_all_non_host_contigs.fasta -c 0.95 -o rr034b1_feces_all_non_host_contigs_DEDUP.fasta -M 5000

```

I ran the above on the UC Berkeley computing cluster. It can be run locally but it is extremeley slow.

---

3. Next,  after the contigs are deduped (or simulataneously as you are doing this), you can download all the (a) nucleotide and (b) protein full genome **reference sequences** under (i) Betacoronavirus (taxid:694002) and (ii) Alphacoronavirus (taxid: 693996) from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). I concatenated the AlphaCoV and BetaCoV nucelotide and protein files such that I had only two CoV reference fasta files: one for nt and one for protein. 

---

4. Now, make sure that the command line version of NCBI-Blast is installed on your computer. See [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for directions. If you are using Mac OS X, running the .dmg installer will probably give you the most success. To test if your installation worked, from anywhere in the command line, try typing 'blastn' and hitting enter. If you receive the following message, then you should be good to go:

```
BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified
```

---

4. Once Blast is installed, you need to make a reference database from your downloaded coronavirus sequences. To do this, on the command line, you can "cd" into the folder where these are contained and use the following command on each of the fasta files produced in step 2 to build these into two different blast databases:

```
makeblastdb –in NCBI_all_CoV_nt.fasta –dbtype nucl –parse_seqids -out CoV_nt
makeblastdb –in NCBI_all_CoV_protein.fasta –dbtype prot –parse_seqids -out CoV_nt
```

Note that you may have to delete and retype the dashes above in your own command line run. This may not copy/paste easily. The commands above should generate a suite of files in the same folder that have the prefix specified after the "out" command.

---

5. Now, that the reference sequence is in hand and the non-host contigs deduplicated, you can kick off a command line blast for each contig subset on the two above databases (so six runs in total). I again ran these on the UC Berkeley computing cluste. Basically, I am re-doing a more focused version of IDseq to see if any other "hits" to CoVs specifically shake out.

