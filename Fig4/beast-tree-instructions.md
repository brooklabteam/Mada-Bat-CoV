## Building a phylogenetic tree using BEAST

### Gathering background sequences from GenBank

The goal of the analysis was to estimate the timing of most recent common ancestor of all Nobecoviruses, as well as estimate the timing of divergence of the Madagascar *Pteropus* clade from the other 4 lineages defined in our paper. To do this, we first compiled all 18 full genome 32 Nobecovirus sequences available on NCBI virus, excepting those two including the p10 orthoreovirus insertion (GCCDC1 lineage). In order to estimate the time of Most Recent Common Ancestor (MRCA) of the Malagasy clades with the corresponding African and Asian clades, we opted to construct a Bayesian timetree.

### Alignment and model selection

After compiling sequences for use in the tree, we renamed the sequences into BEAST format (see script "pre-beast-name.R"), including the accession number and the collection date. For any sequences for which only a collection year was reported, we set the middle of the year (July 15) as the corresponding date. We then aligned these sequences using the [MAFFT](https://mafft.cbrc.jp/alignment/server/) algorithm in Geneious, and evaluated the optimal nucleotide substitution model using [ModelTest-NG](https://github.com/ddarriba/modeltest). This corresponded to GTR+I+G4 for the dataset in hand. 

This folder includes the script used to rename sequences ("pre-beast-name.R"), the original genome sequence file for both the run excluding GCCDC1 ("fullgenome-Nobeco-no-recomb.fasta") and that including it ("allNobeCoVs_BEAST.fasta"), and the multiple sequence alignments (MSA) (in folder "1-genomes-and-alignments") for both. The output from ModelTest-NG for the optimal nucleotide substitution model for this MSA can be accessed in the folder "2-modeltest-ng-output".

---

### Building a phylogenetic tree using BEAST2

The BEAST and BEAST2 community mains a number of extremely helpful tutorials that you should practice to get BEAST up and running on your home computer. See [here](https://taming-the-beast.org/tutorials/). 

After model selection, we next built a Bayesian phylogenetic tree for the *Nobecovirus* phylogeny in BEAST2, using the GTR+I+G4 nucleotide substitution model as recommended by ModelTest-NG. The first step in this process requires generation of a .xml file for input into BEAST, using the program BEAUti. The GTR+I+G4 substitution model is easily specified in BEAUti.  For specifying less common substitution models in BEAUTi, see this great blog post [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/), the recommendations [here](https://groups.google.com/g/ggplot2/c/H50aGubqt2U), or [this rundown](http://www.iqtree.org/doc/Substitution-Models) of comparable model structures for advice. Additionally, it is possible to generate .xml files outside of BEAUti, for example in a pipeline scenario, though this approach was not needed for this phylogeny.

To prepare the .xml file, we used the following parameters in tab inputs at the top of the screen in BEAUti:
 - **Tip Dates**: We used the date of sample collection as the "Tip Date." For any sample from GenBank which only listed year of collection, we set the tip date to July 15 of the year of collection. Each alignment was uploaded to BEAUti with sequence names arranged so as to easily decipher the date.
 - **Site Model**: Following output from ModelTest-NG, we selected a "Gamma Site Model" with Gamma Category 4 and estimated the proportion of invariate sites.
- **Clock Model**: We built comparative phylogenies using  both a strict molecular clock and a relaxed exponential molecular clock, specified within a non-parameteric Bayesian Skyline Coalescent model, following previous approaches for bat coronvirus analyses (i.e. [Lau et al. 2020](https://journals.asm.org/doi/full/10.1128/JVI.02219-09).
- **Priors**: We used a Bayesian Skyline Coalescent model. The clock rate prior was set to a lognormal distribution with a mean of 0.001, following published values for RNA viruses ([Jenkins et al. 2014](https://link.springer.com/article/10.1007/s00239-001-0064-3)), and all other priors were left at default values specified in BEAUti. Both xml files for the strick and relaxed molecular clocks are available in the subfolder "3-BEAST-Nobeco" within the Fig4 folder.
- **MCMC**: We used an MCMC chain length of 700,000,000 iterations and set tracelog and treelog every 10,000 iterations. All other MCMC metrics were left at default. 
- **Population Size**: The Bayesian Skyline Coalescent model by default assumes that the population size of the dataset will change 5 times spaced evenly across the course of this sampling period. Because our available samples were limited and spanned a wide geographic area, we edited this parameter to permit only one population size. You can make this edit in BEAUti by clicking "View" in the top panel, selecting "Show Initialization panel" and then specifying the dimension of "bPopSizes" and "GroupSizes" both to be 1 instead of 5.

All xml files are available for viewing in the folder "3-BEAST-Nobeco".

---

### Visualizing Bayesian timetree

Output from BEAST is available in the sub-folder "BEAST-nobeco-all" and "BEAST-nobeco-recomb. The initial 10% of MCMC iterations were removed as burn-in. Parameter convergence was assessed visually using Tracer v1.6. As has been previously reported (e.g. [Razanjatovo et al. 2015](https://virologyj.biomedcentral.com/articles/10.1186/s12985-015-0271-y)), the Bayesian Skyline tree with a relaxed molecular clock and an exponential distribution offered the best fit to the data. We thus selected this data subset, used TreeAnnotator to average across the BEAST tree output, and visualized the resulting tree in the program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). After checking for basic alignment with parameteric phylogenies generated from RAxML, we then converted the Bayesian tree which is output in Nexus format to Newick format by exporting from FigTree. We then imported the resulting Newick file of the average tree in R and visualized it in the script "Fig4.R".
