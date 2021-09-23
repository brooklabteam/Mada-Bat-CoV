# Mada-Bat-CoV

This github repo contains documentation of analyses for the Madagascar bat coronavirus paper. Goal of this repo is to document analysis to define which samples from fecal and urine datasets have detectable coronavirus signal, to build a phylogeny of those CoVs, and to assemble full genomes for as many as possible. Additionally, we will undertake a few simple ecological analyses to highlight the seasonality of bat virus shedding.

The repo contains sub-directories with scripts and data to make each figure, as well as directions [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/contig-blast-directions.md) for the pipeline used to determine which samples are "positive" for the pathogen of interest (in this case, CoVs).

Within  the Fig. 3 sub-directory, you will find details on how to build a phylogenetic tree [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig3/Phylo-Tree-Directions.md). This script outlines steps used to build Fig 3A and 3B, as well as the amino acid trees in Fig4.

Additionally, you will find details on how to build a Bayesian phylogenetic tree using [BEAST2](http://www.beast2.org/) within the subfolder for 'FigS2' or linked [here](https://github.com/brooklabteam/Mada-Bat-CoV/blob/main/Fig-S2/beast-tree-intstructions.md). We used this analysis to estimate time to most recent common ancestor for *Nobecoviruses* broadly and the newly-defined Madagascar *Pteropus* lineage specifically.

Beyond the details listed above, all scripts to produce each figure should be self-contained and output results to the "final-figures" folder.



