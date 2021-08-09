rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"

setwd(paste0(homewd, "/Fig3"))

#load the fig3a tree
tree <-  read.tree(file = paste0(homewd, "Fig3/B-RdRP-phylogeny/raxml-output/T3.raxml.supportTBE"))


#take a quick look in base R
plot(tree)

#and in ggtree
p1 <- ggtree(tree)
print(p1) #VERY NICE! but it currently does not have all the details shown in color

#root the tree in the outgroup
rooted.tree <- root(tree, which(tree$tip.label == "NC_010800_1_|Turkey_coronavirus__complete_genome"))
ggtree(rooted.tree) #+ geom_tiplab() # rooted by the outgroup C

#check these
tree$node.label
tree$tip.label

p2 <- ggtree(tree) + geom_tiplab() + geom_nodelab()
print(p2)