## Building a phylogenetic tree using BEAST

### Gathering background sequences from GenBank

The goal of the analysis was to estimate the timing of most recent common ancestor of all Nobecoviruses, as well as estimate the timing of divergence of the Madagascar *Pteropus* clade from the other 4 lineages defined in our paper. In the end, we concluded that the depth of sampling over time was too weak to support demographic reconstructions of this kind. However, we here leave instructions for our process to facilitate future work when additional sequences become available.

As it stands, our tree uses a subset of 32 Nobecovirus sequences collected for the maximum likelihood tree visualized in Figure 3B of the main text. In order to estimate the time of Most Recent Common Ancestor (MRCA) of the Malagasy clades with the corresponding African and Asian clades, we opted to additionally construct a Bayesian timetree. See instructions for Figure 3B [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig3/Phylo-Tree-Directions.md) for details on the sequences used in this analysis. Our Bayesian tree corresponded to all but one of the Fig3B sequences which corresponded to the Nobecovirus subgenus (we dropped the one *Eidolon dupreanum* sequence becuase it was so fragmented and compromised the integrity of the rest of the alignment).

### Alignment and model selection

After compiling sequences for use in the tree, we renamed the sequences into BEAST format (see script "prep_BayesianTree.R"), including the accession number and the collection date. For any sequences for which only a collection year was reported, we set the middle of the year (July 15) as the corresponding date. We then aligned these sequences using the [MAFFT](https://mafft.cbrc.jp/alignment/server/) algorithm in Geneious, trimmed them to the length of the longest available sequence (386bp), and evaluated the optimal nucleotide substitution model using [ModelTest-NG](https://github.com/ddarriba/modeltest). This corresponded to TIM2+G4 for the dataset in hand. 

We also tried the same for just the Pteropus clade sequences (only 9), for which ModelTest-NG evaluated TrN+G4 as the best fit model for the data.

This folder includes the script used to rename sequences ("prep_BayesianTree.R"), the original genome sequence file with ("allRdRp-Nobeco-BEAST.fasta"), and the multiple sequence alignments (MSA) (in folder "alignments") both befoe and after renaming to BEAST form and after trimming to even length. The output from ModelTest-NG for the optimal nucleotide substitution model for this MSA (folder "modeltest-ng-output").

---

### Building a phylogenetic tree using BEAST2

The BEAST and BEAST2 community mains a number of extremely helpful tutorials that you should practice to get BEAST up and running on your home computer. See [here](https://taming-the-beast.org/tutorials/). 

After model selection, we next built a Bayesian phylogenetic tree for the Nobecovirus and Pteropus CoV phylogenies in BEAST2, using the TIM2+G4 and TrN+G4 nucleotide substitution models as recommended by ModelTest-NG. The first step in this process requires generation of a .xml file for input into BEAST, using the program BEAUti. The TIM2+G4 substitution model is not easily specified in BEAUti and requires modifications: we represented it by selecting a GTR+G4 model with linked substitution rates for AC with AT and CG with GT, and with CT and AG parameters fixed at 1.  For specifying less common substitution models in BEAUTi, see this great blog post [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/), the recommendations [here](https://groups.google.com/g/ggplot2/c/H50aGubqt2U), or [this rundown](http://www.iqtree.org/doc/Substitution-Models) of comparable model structures. Additionally, it is possible to generate .xml files outside of BEAUti, for example in a pipeline scenario, though this approach was not needed for this phylogeny.

To prepare the .xml file, we used the following parameters in tab inputs at the top of the screen in BEAUti:
 - **Tip Dates**: We used the date of sample collection as the "Tip Date." For any sample from GenBank which only listed year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uploaded to BEAUti with sequence names arranged so as to easily decipher the date.
 - **Site Model**: Following output from ModelTest-NG, we selected a "Gamma Site Model" with Gamma Category 4 and the linked rates described above.
- **Clock Model**: We built comparative phylogenies using  both a strict molecular clock and a relaxed exponential molecular clock, specified within a non-parameteric Bayesian Skyline Coalescent model, following previous approaches for bat coronvirus analyses (i.e. [Lau et al. 2020](https://journals.asm.org/doi/full/10.1128/JVI.02219-09).
- **Priors**: We used a Bayesian Skyline Coalescent model. The clock rate prior was set to a lognormal distribution with a mean of 0.001, following published values for RNA viruses ([Jenkins et al. 2014](https://link.springer.com/article/10.1007/s00239-001-0064-3)), and all other priors were left at default values specified in BEAUti. Both xml files for the strick and relaxed molecular clocks are available in the subfolder "RdRpBEASTxml" within the FigS2 folder.
- **MCMC**: We used an MCMC chain length of 100,000,000 iterations and set tracelog and treelog every 10,000 iterations. All other MCMC metrics were left at default. 
- **Population Size**: The Bayesian Skyline Coalescent model by default assumes that the population size of the dataset will change 5 times spaced evenly across the course of this sampling period. This parameter was left as default in the original runs (the xml files shown here) but then edited in the AllNobecoRelaxed.xml to only permit one population size. In this latter case, you can make this edit in BEAUti by clicking "View" in the top panel, selecting "Show Initialization panel" and then specifying the dimension of "bPopSizes" to be 1 instead of 5.

All xml files are available for viewing in the folder "RdRpBEASTxml".

---

### Visualizing Bayesian timetree

Output from BEAST is available in the sub-folder "beast-out". The initial 10% of MCMC iterations were removed as burn-in. Parameter convergence was assessed visually using Tracer v1.6. As has been previously reported (e.g. [Razanjatovo et al. 2015](https://virologyj.biomedcentral.com/articles/10.1186/s12985-015-0271-y)), the Bayesian Skyline tree with a relaxed molecular clock and an exponential distribution offered the best fit to the data. We thus selected this data subset, used TreeAnnotator to average across the BEAST tree output, and visualized the resulting tree in the program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). The output never converged, as a result of low sample size and a small genome fragment. As such, we ultimately decided to not include this analysis in the final manuscript and did not also plot the output tree in R.

