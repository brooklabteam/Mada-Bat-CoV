# Mada-Bat-CoV

This github repo contains documentation of analyses for the Madagascar bat coronavirus paper. Goal of this repo is to document analysis to define which samples from fecal and urine datasets have detectable coronavirus signal, to build a phylogeny of those CoVs, and to assemble full genomes for as many as possible. Additionally, we will undertake a few simple ecological analyses to highlight the seasonality of bat virus shedding.

---

1. First, we downloaded all non-host contigs derived from fecal, urine, or saliva and mapping to ANY taxon from IDseq in the 'RR034B1_feces' project, the 'RR034B2_urine_wholeblood_novaseq' project, and the 'RR034B_throat_swab_raw_RNASeq_NovaSeq'project. (Note that we did not take contigs from HeLa controls or water, and in the case of the 'RR034B2_urine_wholeblood_novaseq' project, we only looked at urine samples). This can be done in bulk, manually, on IDseq.net in the top right-hand corner.

Note that when you download all the non-host contigs, it will produce a folder with a separate fasta file for each sample, which lists the contigs by node number but does not include the sample ID. Before joining all the contigs (nodes) together, you need to distinguish them by sample ID. I wrote an Rscript that parses this for each filetype (rename-fastas-feces, rename-fastas-urine, rename-fastas-throat). To rename your files and the headers within them, copy the appropriate Rscipt for the tissue type into your appropriate downloads folder, cd into that folder on the command line, and simply type (example here for feces): 

```

Rscript rename-fastas-feces.R

```

Once the file finishes running, concatenate all the abbreviated filenames into a compiled file for downstream analyses. Note that these files are too big to include in the GitHub repo.

---

2. Now, because many of the contigs will overlap in sequences and slow our mapping down, we can deduplicate the compiled .fasta for each of the contigs for all tissues using the program [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/). You'll need to [install this](https://github.com/weizhongli/cdhit/wiki/2.-Installation) on your home computer first. If using MacOS, I recommend using bioconda to do it--see [here](https://anaconda.org/bioconda/cd-hit). Once installed, you can deduplicate each of the contig files via the following command line script in the same folder as the downloaded contigs (here shown just for the fecal contigs file -- you will need to run it for all three tissue types). This command :

```
cd-hit-est -i rr034b1_feces_all_non_host_contigs.fasta -c 0.95 -o rr034b1_feces_all_non_host_contigs_DEDUP.fasta -M 0

```

I ran the above on the UC Berkeley computing cluster (Savio). It can be run locally but it is extremely slow. The fecal sample run took about 5 hours on my personal computer and about 2 hours on Savio. The throat and urine hav fewer contigs and are much faster.

Here's the bash script I used on Savio to run the fecal example:

```
#!/bin/bash
#SBATCH --job-name=cdhit-feces
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

module load gcc/7.4.0
module load cmake/3.15.1

/global/home/users/cbrook/cd-hit/cd-hit-est -i rr034b1_feces_all_non_host_contigs.fasta -c 0.95 -o rr034b1_feces_all_non_host_contigs_DEDUP.fasta -M 0
```

---

3. Next,  after the contigs are deduped (or simulataneously as you are doing this), you can download all the (a) nucleotide and (b) protein full genome **reference sequences** under (i) Betacoronavirus (taxid:694002) and (ii) Alphacoronavirus (taxid: 693996) from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). I concatenated the AlphaCoV and BetaCoV nucelotide and protein files such that I had only two CoV reference fasta files: one for nt and one for protein. 

---

4. Now, you need to make sure that the command line version of NCBI-Blast is installed on your computer. See [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for directions. If you are using Mac OS X, running the .dmg installer will probably give you the most success. To test if your installation worked, from anywhere in the command line, try typing 'blastn' and hitting enter. If you receive the following message, then you should be good to go:

```
BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified
```

Also, note that I ran all the blast searches on Savio, as well. If you go this route, you will need to install Blast to the computing cluster.

---

4. Once Blast is installed, you need to make a reference database from your downloaded coronavirus sequences. To do this, on the command line, you can "cd" into the folder where these are contained and use the following command on each of the fasta files produced in step 2 to build these into two different blast databases:

```
makeblastdb –in NCBI_all_CoV_nt.fasta –dbtype nucl –parse_seqids -out CoV_nt
makeblastdb –in NCBI_all_CoV_protein.fasta –dbtype prot –parse_seqids -out CoV_nt
```

Note that you may have to delete and retype the dashes above in your own command line run. This may not copy/paste easily. The commands above should generate a suite of files in the same folder that have the prefix specified after the "out" command.

---

5. Now that the reference sequence is in hand and the non-host contigs deduplicated, you can kick off a command line blast for each contig subset on the two above databases (so six runs in total). I again ran these on the UC Berkeley computing cluster. Basically, I am re-doing a more focused version of IDseq to see if any other "hits" to CoVs specifically shake out.

You will run six BLASTs in total: 3 nucleotide and 3 protein BLASTs, one for each of the deduplicated contigs above. First, run a "blastn"" alignment of deduplicated set of contigs with the CoV_nt database, then run a "blastx" alignment of the deduplicated set of contigs with the CoV_aa database. Scripts for both are listed below (make sure that you upload ALL the CoV_aa and CoV_nt files into the same folder for this to be able to run):

```
blastn -word_size 10 -evalue 0.001 -query rr034b1_feces_all_non_host_contigs_DEDUP.fasta -db CoV_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out 20210721_Mada_Bat_CoV_blast_nt.txt

blastx -word_size 3 -evalue 0.001 -query rr034b1_feces_all_non_host_contigs_DEDUP.fasta -db CoV_aa -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out 20210721_Mada_Bat_CoV_blast_prot.txt

```

Again, these can be run locally, but I used Savio to speed things up. The searches did not take so long (minutes!) after the deduplication step above for the throat and urine. The fecal search was longer.

Here is the bash script for the deduped fecal contigs for the nucleotide blast:

```
#!/bin/bash
#SBATCH --job-name=blastn-feces
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

module load gcc/7.4.0
module load cmake/3.15.1

/global/home/users/cbrook/ncbi-blast/bin/blastn -word_size 10 -evalue 0.001 -query rr034b1_feces_all_non_host_contigs_DEDUP.fasta -db /global/scratch/cbrook/NCBI-databases/CoV_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out 20210721_Mada_Bat_CoV_blast_feces_nt.txt
```

---

6.  After the blast finishes, you'll want to curate a bit to the high quality hits. After Amy's lead, I went ahead and parsed for alignments that show alignment length > 100 aa and bit score > 100.

Here's the script for the nt parse (feces as example):

```
cat 20210721_Mada_Bat_CoV_blast_feces_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > 20210721_Mada_Bat_CoV_blast_feces_nt_results_100len5eval.txt


```

And for the protein parse (feces again):

```

cat 20210721_Mada_Bat_CoV_blast_feces_prot.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > 20210721_Mada_Bat_CoV_blast_feces_prot_results_100len100bit.txt


```

---

7. Once the blast results have beeen sub-selected a bit, you can summarize them to link back the hits to the samples of interest. Within the same folder as your output, try the following script to save the unique contig IDs which align to CoVs (example here for blastn alignment of fecal samples):

```
cat 20210721_Mada_Bat_CoV_blast_feces_nt_results_100len5eval.txt | awk '{print $1}' | sort | uniq > 20210721_Mada_Bat_CoV_unique_contigs_feces_nt.txt

```

And here to save the unique sample IDs for the same example:

```
cat 20210721_Mada_Bat_CoV_unique_contigs_feces_nt.txt | awk -F\_ '{print $1"_"$2}' | sort | uniq > 20210721_Mada_Bat_CoV_unique_sampleID_feces_nt.txt

```

---

8. And do the same for the other sample types and for the blastx outputs. Load the outputs into R and determine which samples with meta-data were infected at various times/places. I put a bash script ('curate-blast-output.txt') in the 'blast-output' folder that runs through steps 6 through 8 for all of the blast output from this project. You can run it with:

```
sh -e curate-blast-output.txt
```

---

9. Now, import the curated contigs into R and compare them against the IDseq hits for what is CoV positive and how it maps into the meta-data. It looks like no throat samples were CoV hits, as is consistent with what is recovered from IDseq. See R script for further comparison of manual pipeline hits for fecal and urine samples.

