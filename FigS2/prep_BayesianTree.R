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
dat$collection_date <- as.Date(dat$collection_date)
dat$beast_name <- paste0(dat$accession_num, "-", dat$collection_date)

#and load the aligned -first all Nobeco

RdRp <- read.alignment(file = paste0(homewd, "FigS2/alignments/AlignAllNobecoRdRp_Extract.fasta"), format="fasta", forceToLower=FALSE)
head(RdRp)

#check
setdiff(RdRp$nam, dat$tip_label)
setdiff(dat$tip_label, RdRp$nam) #many

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
write.fasta(RdRp$seq, names = dat.new$beast_name, file.out = paste0(homewd, "FigS2/alignments/AlignAllNobecoRdRp_Extract_BEAST.fasta"), open="w", as.string = T)

#then send to modeltest, BEAUti, and BEAST

#and the Pteropus - cut
# 
# Pter<- read.alignment(file = paste0(homewd, "FigS2/alignments/Align_Madagascar_Pteropus_Subclade_RdRp_Extract.fasta"), format="fasta", forceToLower=FALSE)
# head(Pter)
# 
# 
# #check
# setdiff(Pter$nam, dat$tip_label)
# setdiff(dat$tip_label, Pter$nam) #many
# 
# #and rename
# dat.orig2 <- cbind.data.frame(tip_label = Pter$nam)
# head(dat.orig2)
# tail(dat.orig2)
# 
# #and merge
# dat.new2 <- merge(dat.orig2, dat, by = "tip_label", all.x = T, sort = F)
# head(dat.new2)
# tail(dat.new2)
# 
# 
# Pter$seq  <- lapply(Pter$seq, toupper)
# 
# #now resave the files with beast name
# write.fasta(Pter$seq, names = dat.new2$beast_name, file.out = paste0(homewd, "FigS2/alignments/Align_Madagascar_Pteropus_Subclade_RdRp_Extract_BEAST.fasta"), open="w", as.string = T)
# 
# 
# #then send to modeltest, BEAUti, and BEAST
# 
# 
