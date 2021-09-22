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

dat <- read.csv(file="figS2_Bayesian_tree_RdRp_metadata.csv", header = T, stringsAsFactors = F)
head(dat)
dat$collection_date <- as.Date(dat$collection_date, format = "%Y-%m-%d")
dat$beast_name <- paste0(dat$accession_num, "-", dat$collection_date)

#and load the aligned 
RdRp <- read.alignment(file = paste0(homewd, "FigS2/AlignRdRP-Full-NoGamma-RdRp.fasta"), format="fasta", forceToLower=FALSE)
head(RdRp)

#check
setdiff(RdRp$nam, dat$tip_label)
setdiff(dat$tip_label, RdRp$nam)

#and rename
dat.orig <- cbind.data.frame(tip_label = RdRp$nam)
head(dat.orig)
tail(dat.orig)

#and merge
dat.new <- merge(dat.orig, dat, by = "tip_label", all.x = T, sort = F)
head(dat.new)
tail(dat.new)


RdRp$seq  <- lapply(RdRp$seq, toupper)

#now resave the files with beast name
write.fasta(RdRp$seq, names = dat.new$beast_name, file.out = "AlignRdRP-Full-NoGamma-RdRp-BEAST.fasta", open="w", as.string = T)

#then send to modeltest, BEAUti, and BEAST

