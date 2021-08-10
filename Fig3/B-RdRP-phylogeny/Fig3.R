rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"

setwd(paste0(homewd, "/Fig3"))

#load the fig3a tree
treeA <-  read.tree(file = paste0(homewd, "Fig3/A-full-genome-phylogeny/Fig3A-raxml-output/Fig3A.raxml.supportTBE"))


#take a quick look in base R
plot(treeA)

#and in ggtree
p1 <- ggtree(treeA)
print(p1) #VERY NICE! but it currently does not have all the details shown in color

#root the tree in the outgroup
rooted.tree.A <- root(treeA, which(tree$tip.label == "NC_010800_1_|Turkey_coronavirus__complete_genome"))
ggtree(rooted.tree.A) #+ geom_tiplab() # rooted by the outgroup C

#check these
treeA$node.label
treeA$tip.label

p2 <- ggtree(treeA) + geom_tiplab() + geom_nodelab()
print(p2)


#### And Fig 3 B
