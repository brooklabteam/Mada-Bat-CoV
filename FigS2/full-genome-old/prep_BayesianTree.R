rm(list=ls())
#load the tree data
library(plyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(seqinr)


#load trees and make into multipanel amino acid tree plot

homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"
setwd(paste0(homewd, "FigS2")) 

dat <- read.csv(file="figS2_Bayesian_tree_metadata.csv", header = T, stringsAsFactors = F)
head(dat)

dat$beast_name <- paste0(dat$accession_num, "-", dat$collection_date)

#and load the sequences
all.beta <- read.fasta(file = paste0(homewd, "FigS2/allBetaCoVsBEAST.fasta"),forceDNAtolower = F, as.string = T)
head(all.beta)

#check
setdiff(attr(all.beta, "name"), dat$tip_label)
setdiff(dat$tip_label, attr(all.beta, "name"))

#and rename
dat.orig <- cbind.data.frame(tip_label = attr(all.beta, "name"))
head(dat.orig)
tail(dat.orig)

#and merge
dat.new <- merge(dat.orig, dat, by = "tip_label", all.x = T, sort = F)
head(dat.new)
tail(dat.new)


#all.beta  <- lapply(all.beta$seq, toupper)

#now resave the files with beast name
write.fasta(all.beta, names = dat.new$beast_name, file.out = "AllBetaCoV_BEASTnames.fasta", open="w", as.string = T)

#now align in MAFFT, send to BEAUti, and send to BEAST

#and subset just Nobecoviruses
rm(list=ls())

homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"
setwd(paste0(homewd, "FigS2")) 

dat <- read.csv(file="figS2_Bayesian_tree_metadata.csv", header = T, stringsAsFactors = F)
head(dat)

dat$beast_name <- paste0(dat$accession_num, "-", dat$collection_date)

#and load the sequences
all.beta <- read.fasta(file = paste0(homewd, "FigS2/allBetaCoVsBEAST.fasta"),forceDNAtolower = F, as.string = T)
head(all.beta)

#check
setdiff(attr(all.beta, "name"), dat$tip_label)
setdiff(dat$tip_label, attr(all.beta, "name"))

#and rename
dat.orig <- cbind.data.frame(tip_label = attr(all.beta, "name"))
head(dat.orig)
tail(dat.orig)

#and merge
dat.new <- merge(dat.orig, dat, by = "tip_label", all.x = T, sort = F)
head(dat.new)
tail(dat.new)

dat.nobec = subset(dat.new, sub_group=="Nobecovirus") #17
list.nobec <- rownames(dat.nobec)
nobec.seq <- list()

for (i in 1:length(list.nobec)){
  id = as.numeric(list.nobec[i])
  tmp.seq = all.beta[[id]]
  nobec.seq[[i]] <- tmp.seq
}

#and save with beast name
write.fasta(nobec.seq, names = dat.nobec$beast_name, file.out = "NobeCoV_BEAST.fasta", open="w", as.string = T)

