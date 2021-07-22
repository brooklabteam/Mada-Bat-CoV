rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

#now load the hits for CoVs (urine and fecal only -- throat had 0 hits) and summarise

homewd = getwd()

dat.fec.prot <- read.delim(file = paste0(homewd,"/blast-output/20210721_Mada_Bat_CoV_unique_sampleID_feces_prot.txt"), header = F)
dat.fec.nt <- read.delim(file = paste0(homewd,"/blast-output/20210721_Mada_Bat_CoV_unique_sampleID_feces_nt.txt"), header = F)
dat.ur.prot <- read.delim(file = paste0(homewd,"/blast-output/20210721_Mada_Bat_CoV_unique_sampleID_urine_prot.txt"), header = F)
dat.ur.nt <- read.delim(file = paste0(homewd,"/blast-output/20210721_Mada_Bat_CoV_unique_sampleID_urine_nt.txt"), header = F)

#and summarise
dat.ur.prot <- data.frame(sampleID=dat.ur.prot[,1], tissue = "urine", blast = "protein")
dat.ur.nt <- data.frame(sampleID=dat.ur.nt[,1], tissue = "urine", blast = "nt")
dat.fec.prot <- data.frame(sampleID=dat.fec.prot[,1], tissue = "feces", blast = "protein")
dat.fec.nt <- data.frame(sampleID=dat.fec.nt[,1], tissue = "feces", blast = "nt")

all.hits <- rbind(dat.ur.prot, dat.ur.nt, dat.fec.prot, dat.fec.nt)

#all unique is a vector of unique sampleIDs which are CoV (+) by this pipeline
all.unique <- c(ddply(all.hits, .(sampleID), summarise))

#now load our original dataset and compare those that match as positive.

