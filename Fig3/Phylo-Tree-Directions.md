

## Building a phylogenetic tree

This tutorial outlines methods for building a maximum-likelihood phylogenetic tree using the program RAxML. We will build (A) a full genome phylogeny of all full genome sequences of Alpha- and Betacoronaviruses on [NCBI-Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), fitting our new Madagascar whole genome sequences into the picture, and (B) a partial tree focused on just the RdRp gene of the virus. I outline directions below for how to build (A) first.

---

### A. Full-genome Alpha- BetaCoV phylogeny

1. I first went to NCBI-virus (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and selected all **full genome** sequences available for Betacoronavirus (tax id 694002), which brought up 412,467 hits. I did the same for Alphacoronavirus (tax id 693996), which brought up 1,246 hits. For each CoV genus, I downloaded a .csv file which I processed in the attached R script to select viruses which met the following criteria:
 - all full genome betaCoVs with a bat host (98 in total)
 - all full genome betaCoV **reference** sequences with a non-bat host (14)
 - all full genome alphaCoVs with a bat host (118 in total)
 - all full genome alphaCoV **reference** sequences with a non-bat host (13)
 
 This led to a total of 243 records. I downloaded these into a folder on my home computer in bulk from NCBI using the following code loaded into my web browser (produced from the Rscript):
 
 
```
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=KY417142, MN611520, MG596802, MG596803, KJ473820, KJ473812, KJ473811, KJ473813, MK211374, MK211375, MK211376, MK211377, MK211378, KJ473815, KJ473814, KJ473816, MK211379, KJ473822, JX993988, KU973692, MT350598, NC_030886, KU762338, MN611519, KF569996, KX442565, MN996532, DQ412042, KY417145, MW719567, DQ412043, MG762674, JX993987, DQ071615, KC881006, KY417143, KY417144, KY417146, KY417147, KY417148, KY417149, KY417150, FJ588686, KY417151, KY417152, KC881005, MW218395, KF367457, KT444582, KX442564, KP886808, KP886809, NC_025217, KF636752, MG772933, MG772934, NC_014470, NC_009019, NC_009020, NC_009021, MZ081376, MZ081377, MZ081378, MZ081379, MZ081380, MZ081381, MZ081382, MZ190137, MZ190138, LC556375, KY352407, MF593268, KJ473821, KC869678, HM211098, HM211099, HM211100, HM211101, GU190215, EF065505, EF065506, EF065507, EF065508, EF065509, EF065510, EF065511, EF065512, EF065513, EF065514, EF065515, EF065516, MH002340, MH002341, MH002337, MH002338, MH002339, MH002342, DQ648794, NC_003045, NC_039207, NC_038294, NC_019843, NC_006577, NC_004718, NC_045512, NC_048217, NC_006213, NC_026011, AC_000192, NC_017083, NC_012936, NC_001846, MN611522, MN611523, MN611524, JQ989271, NC_018871, JQ989270, MN611517, MN611518, MT663548, KJ473795, KJ473799, KJ473797, KJ473800, KJ473798, KJ473796, NC_028811, KJ473806, KJ473810, NC_028833, KJ473809, NC_028814, KJ473807, NC_028824, KJ473808, MK211373, MK211369, MK211370, MK211371, MK211372, MN611525, MN611521, JQ989269, JQ989267, JQ989268, JQ989272, JQ989273, JQ989266, NC_048216, NC_046964, NC_032107, NC_022103, MZ081383, MZ081384, MZ081385, MZ081386, MZ081387, MZ081388, MZ081389, MZ081390, MZ081391, MZ081392, MZ081393, MZ081394, MZ081395, MZ081396, MZ081397, MZ081398, MZ081399, MN065811, MK720944, MK720945, MK720946, MH938448, MH938449, MH938450, MH687934, MH687935, MH687936, MH687937, MH687938, MH687939, MH687940, MH687941, MH687942, MH687943, MH687944, MH687945, MH687946, MH687947, MH687948, MH687949, MH687950, MH687951, MH687952, MH687953, MH687954, MH687955, MH687956, MH687957, MH687958, MH687959, MH687960, MH687961, MH687962, MH687963, MH687964, MH687965, MH687966, MH687967, KY799179, KY073744, KY073745, KY073746, KY073747, KY073748, KF430219, NC_010437, NC_010438, NC_009988, NC_009657, MG742313, EU420137, EU420138, EU420139, EF203064, EF203065, EF203066, EF203067, NC_002306, NC_034972, NC_030292, NC_032730, NC_038861, NC_035191, NC_048211, NC_028752, NC_028806, NC_023760, NC_005831, NC_003436, NC_002645
```

2. Next, I downloaded one I added our three full-genome Madagascar bat-CoVs from IDseq (Project RR034B1_feces). These are the contig names:
 - RR034B_010_NODE_1_length_29122_cov_127.557574 (*Pteropus rufus*, AMB130)
 - RR034B_232_NODE_1_length_28926_cov_14.641786 (*Rousettus madagascariensis*, )
 - RR034B_288_NODE_2_length_28980_cov_27.781642  (*Rousettus madagascariensis*)
 
3. Finally, I wanted to include a full-genome CoV that was neither an AlphaCoV nor a BetaCoV as an outgroup in my RAxML tree, so I picked a bird GammaCoV from NCBI (Turkey coronavirus, Accession Number: NC_010800.1). I placed these all in the same folder as my Alpha/BetaCoV references above, then "cd"ed into the folder on the command line and concatenated with:

```
cat *fasta > AlphaBetaFullGenomePhylo_CoVs.fasta
```

This .fasta file is included in our GitHub repo.

4. From there, I uploaded the concentenated fasta into the [MAFFT program online](https://mafft.cbrc.jp/alignment/server/) for alignment, using default parameter values. Note that the alignment will take around 30-60 minutes in this case. After the alignment returns, I downloaded it as a .fasta file and saved to my computer. I also saved my version of the alignment file ("multiple sequence alignment," or MSA) in the GitHub repo. I then used R (see "prepre_Fig3a.R") to edit the names of each sequence in the MSA, since RAxML won't accept spaces, periods, dashes, slashes, colon, semicolons, or parentheses in the sequence names. 

5. Next, I used the MSA to compare nucleotide substitution models in the program [ModelTest-NG](https://github.com/ddarriba/modeltest). I ran this on the command line. See GitHub documentation for how to set it up on your computing cluster or personal computer (Note: these phylogenies are sufficiently large that running it on your desktop/laptop is not advisable). After uploading this file to the Berkeley computing cluster, Savio, I kicked off ModelTest-NG using the following batch script:

```
#!/bin/bash
#SBATCH --job-name=batCoV
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00


module load vim/7.4 
module load java/1.8.0_121
module load emacs/25.1 
module load cmake/3.15.1
module load python/3.6 
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

export CC=`which gcc`
export CXX=`which c++`


./modeltest-ng -i batCoV/alignment_fullgenomeCoVs_8_7.fasta -t ml -p 8
```

Note that if you want to **test** this process on your home computer, you can do so by simply reducing the number of genomes in the alignment. I recommend 20 or fewer in this case!

6. Once ModelTest-NG finishes (it will take several hours), it is time to build a maximum likelihood tree using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/). See documentation on their website for how to get this running on your home computer and/or computing cluster. In my case, ModelTest-NG told me that the best support was recovered for a "GTR+I+G4" nt substitution model, so that is what I ran in raxml. I followed the [RAxML tutorial](https://github.com/amkozlov/raxml-ng/wiki/Tutorial) to build a tree from my MSA, first checking that RAxML could read the alignment:

```
/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --check --msa alignment_fullgenomeCoVs_8_7_rename.fasta --model GTR+I+G4 --prefix T1
```

...then parsing the alignment to find the appropriate number of threads (12) with which to run RAxML:

```
/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --parse --msa alignment_fullgenomeCoVs_8_7_rename.fasta --model GTR+I+G4 --prefix T2
```

Finally, I kicked off RAxML (including bootstraps) with the following script:

```
#!/bin/bash
#SBATCH --job-name=batcov
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --all --msa alignment_fullgenomeCoVs_8_7_rename.fasta --model GTR+I+G4 --prefix T3  --seed 12 --threads 12 --bs-metric fbp,tbe
```

Once RAxML finished (>12 hours later), I imported the resulting tree into R and built a phylogenetic tree for Fig3A. Note that for quick viewing, you can easily view the tree in the opensource program, [FigTree](http://tree.bio.ed.ac.uk/software/figtree/).

---

### B. Partial genome RdRp from positive samples

To make Fig3B, a partial RdRP phylogeny from all 30 positive samples, I did the following:

1. Downloaded the two reference Nobecovirus full genomes from GenBank (accession numbers NC_030886.1 and NC_009021.1), as well as all the Madagascar RdRp fragments cited in [Razanajatovo et al. 2015](https://doi.org/10.1186/s12985-015-0271-y). Here is the script to download from a web browser:

```
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_030886.1, NC_009021.1, KP696745, KP696741, KP696746, KP696742, KP696743, KP696744, KP696747, KF859764, KF859759, KF859768, KF859771, KF859761, KF859767, KF859763, KF859760, KF859762, KF859758,  KF859770, KF859765, KF859766, KF859769
```

I then joined these with the three full genomes of the Madagascar bats in a single fasta ("Mada_RdRp_All.fasta") and loaded and aligned them on Geneious.

2. I next wanted to map as many contigs from positive samples to these genomes as possible, to identify which of them mapped to the RdRp region. We had way too many contigs to start with, however, so I first blasted contigs from all positive samples to an RdRp database to determine candidates.

3. To accomplish the blast, I first went back to the de-duped contigs from the homepage (results from CD-HIT) and sub-sampled them to **include just contigs from those 30 bats positive for CoV infection**. See script "preprep_Fig3B.R" in this folder for details.

4. Then, I subsampled the alignment described in #1 above to just the RdRp fragment of all genomes and used this to make an nt blast database for the bat CoV RdRp gene. I also translated the subsampled fragments in Geneious and exported them tomake an equivalent file relevant to the protein blast. Once I had the desired fragments into files ("RdRp-extraction.fasta" and "RdRp-extraction-translation.fasta"), I made these into, respectively, a blastn nt and blastx aa database using the following script:

```
makeblastdb –in RdRp-extraction.fasta –dbtype nucl –parse_seqids -out CoV_RdRP_nt
```

And the protein database:

```
makeblastdb –in RdRp-extraction-translation.fasta –dbtype prot –parse_seqids -out CoV_RdRP_aa
```

Note that you may have to retype the commands above into your terminal, as they are very sensitive to copy/paste.

5. Once I had the above RdRp database, I blasted all the non-host contigs from the positives against this, using the following script (I ran this locally and it was done in seconds):

```
blastn -word_size 10 -evalue 0.001 -query all_CoV_pos_contigs_feces.fasta -db /Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/Fig3/B-RdRP-phylogeny/sequences/RdRp-database/CoV_RdRP_nt -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle'  -max_target_seqs 10 -out 20210808_Mada_Bat_CoV_RdRp_blast_feces_nt.txt
```

And the protein blast too:

```
blastx -word_size 3 -evalue 0.001 -query all_CoV_pos_contigs_feces.fasta -db /Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/Fig3/B-RdRP-phylogeny/sequences/RdRp-database-protein/CoV_RdRp_aa -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out 20210808_Mada_Bat_CoV_RdRp_blast_feces_prot.txt
```

I then subselected for high quality hits using this, for nt:


```
cat 20210808_Mada_Bat_CoV_RdRp_blast_feces_nt.txt | awk -F\t '($4>99 && $5<0.00001) {print $1,$3, $4, $5, $8,$9}' > 20210808_Mada_Bat_CoV_RdRp_blast_feces_nt_results_100len5eval.txt

```

and for protein:


```
cat 20210808_Mada_Bat_CoV_RdRp_blast_feces_nt_results_100len5eval.txt | awk -F\t '($4>99 && $6>99) {print $1,$3, $4,$5,$8,$9}' > 20210808_Mada_Bat_CoV_RdRp_blast_feces_prot_results_100len5eval.txt
```

And then I summarized the above into a file containing unique contig ids using this script, for feces:

```
cat 20210808_Mada_Bat_CoV_RdRp_blast_feces_nt.txt | awk '{print $1}' | sort | uniq > 20210808_Mada_Bat_CoV_RdRp_blast_feces_unique_contigs_hiqual.txt

```

And this for protein:

```
cat 20210808_Mada_Bat_CoV_RdRp_blast_feces_prot_results_100len5eval.txt | awk '{print $1}' | sort | uniq > 20210808_Mada_Bat_CoV_RdRp_unique_contigs_feces_prot_hiqual.txt

```

6. With the above list, I then returned to script "preprep_Fig3.R" and collected contigs which were positive hits. Using those positive contigs, I then did an MSA linking each contig in turn to its closest full genome neighbor in Geneious, then selected those which mapped to an overlapping region of the RdRp gene. I also manually went through all contigs produced from the IDseq pipeline which mapped to any coronavirus and did the same mapping exercise. All-in-all, I found seven genome segments (2x P. rufus, 4 x R. madagascariensis, 1x E. dupreanum) which included the RdRp gene, which I included in the Fig3B phylogeny (see file 'all_RR034B_samples_RdRp.fasta').

7. After selecting the RdRp genomes from the Madagascar project, I combined them with (A) all the Madagascar bat RdRp genome fragments that overlapped from [Razanajatovo et al. 2015](https://doi.org/10.1186/s12985-015-0271-y) (these were the 993bp fragments starting with "KP"), (b) all the Madagascar Betacoronavirus sequences from [Joffrin et al. 2019](https://doi.org/10.1038/s41598-020-63799-7) (these I ended up deleting because they did not overlap), (c) all Nobecovirus RdRp fragments cited in Joffrin et al. 2019, and (d) all Betacoronavirus full genome reference genomes (including two Nobecoviruses, a single Hibecovirus, two Sarbecoviruses, one Merbecovirus, and two Embecoviruses).

Here is the script I used to download the Nobecovirus RdRp fragments:

```
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=KX285095,GU065433,HQ728482,KX285299,KX285008,KX285087,KX520659,KX284906,KU182986,KX284911,NC009022,AB683971,AB539082,KU182962
```

8. I aligned all sequences from (7) in Geneious using MAFFT, then sub-selected those that actually overlapped for a good region of RdRp. The final set of 36 sequences used in my RdRp phylogeny can be found in the "all_RdRp_Bat_CoV.fasta" file in the "final-set" subfolder of the "final-RdRp-alignment". Once these were selected, I aligned them again ("Align-Full-Seq-RdRp-8-10-Final.fasta"), then trimmed the alignment to the longest sequence. This cut the alignment to a meager 259 bp ("Align-Severe-RdRp-8-10-Extraction.fasta").

9. I then queried the best nt substitution model on the severe alignment using ModelTest-NG (it was 'TVM+I+G4'), and built a bootstraped maximum likelihood phylogeny using RAxML. I followed the instructions above in part A, #6 to check and parse the alignment before kicking it off. Note that both ModelTest and RAxML went MUCH faster on these genomes.

10. I imported the resulting phylogeny into R to make the tree seen in Figure 3B (see "Fig3.R" script).
