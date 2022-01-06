rm(list=ls())

library(plyr)
library(dplyr)
library(seqinr)


#set wd
homewd = "/Users/carabrook/Developer/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig4/"))

#load the dataset and query
dat <- read.csv(file = "fig4beast_nobecov_metadata_manual.csv", header = T, stringsAsFactors = F)

head(dat)

#now, get the fasta file 
fasta.dat <- read.fasta(file= paste0(homewd,"/Fig4/1-genomes-and-alignments/fullgenome-Nobeco.fasta"), forceDNAtolower = F, as.string = T)

names(fasta.dat)

fasta.meta <-  cbind.data.frame(tip_label = names(fasta.dat))

#add beast name to the main dataset
dat$collection_date <- as.Date(dat$collection_date)
#dat$collection_date <- as.character(dat$collection_date)
#dat$collection_date[is.na(dat$collection_date)] <- paste0(dat$collection_year[is.na(dat$collection_date)], "07-31")
#dat$collection_date <- as.Date(dat$collection_date)

dat$beast_name <- paste0(dat$accession_num, "_", dat$collection_date)

dat.merge <- dplyr::select(dat, tip_label, beast_name)


setdiff(dat.merge$tip_label, fasta.meta$tip_label)
setdiff( fasta.meta$tip_label,dat.merge$tip_label)

#fasta.meta$tip_label[fasta.meta$tip_label=="DQ648794_1_Bat_coronavirus_(BtCoV_133_2005)"] <- "DQ648794_1_Bat_coronavirus__BtCoV_133_2005"

#and add to data
fasta.meta <- merge(fasta.meta, dat.merge, all.x = T, by="tip_label", sort = F)
head(fasta.meta)
fasta.meta$beast_name
#and write over:
write.fasta(fasta.dat, names=fasta.meta$beast_name, file =paste0(homewd,"/Fig4-new/1-genomes-and-alignments/allNobeCoVs_BEAST.fasta") )


#and the version with no recombination
rm(list=ls())


#set wd
homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig4-new/"))


#load the dataset and query
dat <- read.csv(file = "fig4beast_nobecov_metadata_manual.csv", header = T, stringsAsFactors = F)

head(dat)

#now, get the fasta file 
fasta.dat <- read.fasta(file= paste0(homewd,"/Fig4-new/1-genomes-and-alignments/fullgenome-Nobeco-no-recomb.fasta"), forceDNAtolower = F, as.string = T)

names(fasta.dat)

fasta.meta <-  cbind.data.frame(tip_label = names(fasta.dat))

#add beast name to the main dataset
dat$collection_date <- as.Date(dat$collection_date)
#dat$collection_date <- as.character(dat$collection_date)
#dat$collection_date[is.na(dat$collection_date)] <- paste0(dat$collection_year[is.na(dat$collection_date)], "07-31")
#dat$collection_date <- as.Date(dat$collection_date)

dat$beast_name <- paste0(dat$accession_num, "_", dat$collection_date)

dat.merge <- dplyr::select(dat, tip_label, beast_name)


setdiff(dat.merge$tip_label, fasta.meta$tip_label)
setdiff( fasta.meta$tip_label,dat.merge$tip_label)

#fasta.meta$tip_label[fasta.meta$tip_label=="DQ648794_1_Bat_coronavirus_(BtCoV_133_2005)"] <- "DQ648794_1_Bat_coronavirus__BtCoV_133_2005"

#and add to data
fasta.meta <- merge(fasta.meta, dat.merge, all.x = T, by="tip_label", sort = F)
head(fasta.meta)
fasta.meta$beast_name
#and write over:
write.fasta(fasta.dat, names=fasta.meta$beast_name, file =paste0(homewd,"/Fig4-new/1-genomes-and-alignments/NobeCoVs_no_recomb_BEAST.fasta") )

