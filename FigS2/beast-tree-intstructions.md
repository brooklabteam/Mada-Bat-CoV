## Building a phylogenetic tree using BEAST

### Gathering background sequences from GenBank

This tree uses the same subset of sequences collected for the maximum likelihood tree visualized in Figure 3B of the main text. In order to estimate the time of Most Recent Common Ancestor (MRCA) of the Malagasy clades with the corresponding African and Asian clades, we opted to additionally construct a Bayesian timetree. See instructions for Figure 3B [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig3/Phylo-Tree-Directions.md) for details on the sequences used in this analysis. Note that we left off the Turkey Gammacoronavirus outgroup in this case, as it is not necessary to include in a Bayesian timetree.

### Alignment and model selection

After compiling sequences for use in the tree, we renamed the sequences in BEAST format (see script "prep_BayesianTree.R"), including the accession number and the collection date. For any sequences for which only a collection year was reported, we set the middle of the year (July 15) as the corresponding date. We then aligned these sequences using the [MAFFT](https://mafft.cbrc.jp/alignment/server/) algorithm in Geneious, and evaluated the optimal nucleotide substitution model using [ModelTest-NG](https://github.com/ddarriba/modeltest). This corresponded to TPM2uf+I+G4 for the dataset in hand. 

This folder includes the script used to rename sequences ("prep_BayesianTree.R"), the original genome sequence file with BEAST-structured sequence names ("allRdRP-all-Nobeco-BEAST.fasta"), the multiple sequence alignment (MSA) of these genomes prior to trimming ("AlignRdRP-Full-NoGamma-RdRp.fasta"), the trimmed MSA to the 259 bp RdRp fragment ("AlignRdRP-Full-NoGamma-RdRp-BEAST.fasta"), and the output from ModelTest-NG for the optimal nucleotide substitution model for this MSA (folder "modeltest-ng-output").

---

### Building a phylogenetic tree using BEAST2

The BEAST and BEAST2 community mains a number of extremely helpful tutorials that you should practice to get BEAST up and running on your home computer. See [here](https://taming-the-beast.org/tutorials/). 

After model selection, we next built a Bayesian phylogenetic tree for the RdRp CoV phylogeny in BEAST2, using the TPM2uf+I+G4 nucleotide substitution model recommended by ModelTest-NG. The first step in this process requires generation of a .xml file for input into BEAST, using the program BEAUti. The TPM2uf+I+G4 substitution model is not easily specified in BEAUti and requires modifications: we represented it by selecting a GTR+I+G4 model with linked substitution rates for AC with AT, CG with GT, and AG with CT.  For specifying less common substitution models in BEAUTi, see this great blog post [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/) or the recommendations [here](https://groups.google.com/g/ggplot2/c/H50aGubqt2U). Additionally, it is possible to generate .xml files outside of BEAUti, for example in a pipeline scenario, though this approach was not needed for this phylogeny.

To prepare the .xml file, we used the following parameters in tab inputs at the top of the screen in BEAUti:
 - **Tip Dates**: We used the date of sample collection as the "Tip Date." For any sample from GenBank which only listed year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uploaded to BEAUti with sequence names arranged so as to easily decipher the date.
 - **Site Model**: Following output from ModelTest-NG, we selected a "Gamma Site Model" with Gamma Category 4 and the linked rates described above.
- **Clock Model**: We built comparative phylogenies using  both a strict molecular clock and a relaxed exponential molecular clock, specified within a non-parameteric Bayesian Skyline Coalescent model, following previous approaches for bat coronvirus analyses (i.e. [Lau et al. 2020](https://journals.asm.org/doi/full/10.1128/JVI.02219-09).
- **Priors**: We used a Bayesian Skyline Coalescent model. The clock rate prior was set to a lognormal distribution with a mean of 0.001, following published values for RNA viruses ([Jenkins et al. 2014](https://link.springer.com/article/10.1007/s00239-001-0064-3)), and all other priors were left at default values specified in BEAUti. Both xml files for the strick and relaxed molecular clocks are available in the subfolder "RdRpBEASTxml" within the FigS2 folder.
- **MCMC**: We used an MCMC chain length of 100,000,000 iterations and set tracelog and treelog every 10,000 iterations. All other MCMC metrics were left at default. 

---

### Visualizing Bayesian timetree

Output from BEAST is available in the sub-folder "beast-out". The initial 10% of MCMC iterations were removed as burn-in. Parameter convergence was assessed visually using Tracer v1.6. Once verified, we used TreeAnnotator to average across the BEAST tree output, then visualized the resulting tree in R. The script for preparation of the FigS2 timetree can be found as "FigS2.R" in this folder. 

