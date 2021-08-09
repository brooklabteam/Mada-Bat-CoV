rm(list=ls())
## purpose of script is to sub-select contigs of positives to 
## then blast to the RdRp reference database

library(plyr)
library(dplyr)
library(seqinr)

#setwd
homewd <- "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig3/B-RdRp-phylogeny"))

#load all positive fecal contigs
fec_contigs <- read.fasta( file = paste0(homewd, "/Non-Host-Contigs/rr034b1_feces_all_non_host_contigs_DEDUP.fasta"), as.string = T, forceDNAtolower = F)

#get the names
names_fec <- names(fec_contigs)

#load all positive bats and select those contigs which belong to them
dat <- read.csv(file = paste0(homewd, "/metadata/all_NGS_8_3_2021_distribute.csv"), header = T, stringsAsFactors = F)
#get the positives
pos.dat <- subset(dat, CoV==1)
head(pos.dat)
#get the samply names of the positves
names_pos <- paste(sapply(strsplit(pos.dat$czb_id, "_"), '[',1), sapply(strsplit(pos.dat$czb_id, "_"), '[',2), sep = "_")

#select any of the contigs that begin with any of these
#write a function to do this
select.contigs <- function(contig_name, all_pos){
  
  #split the contig name and drop to just the sample name
  sample_name <- sapply(strsplit(contig_name, "_NODE"), '[',1)
  
  #check if included in all_pos
  if(length(intersect(sample_name,all_pos))>0){
    #if yes, return the contig
    return(contig_name)
  }
  #if not, do nothing
}

contigs_from_pos <- lapply(as.list(names_fec), select.contigs,  all_pos=names_pos)

contigs_from_pos <- c(unlist(contigs_from_pos))

#now select those sequences from the list
select.seq <- function(seq1, all_pos_contigs){
 
  name_contig = attr(seq1, "name")
  
  #now select that sequence if it matches
  if(length(intersect(name_contig, all_pos_contigs))>0){
    #and return it
    return(seq1)
  }
  
  #otherwise, leave off

}

positive_contigs <- lapply(fec_contigs,select.seq, all_pos_contigs=contigs_from_pos)

#remove NULL entries
positive_contigs <- positive_contigs[lengths(positive_contigs) > 0]

#and save
names_pos = names(positive_contigs)
write.fasta(positive_contigs, names = names_pos, file.out = paste0(homewd, "/Fig3/B-RdRp-phylogeny/blast-pos-contigs/all_CoV_pos_contigs_feces.fasta"), as.string = T, open="w")


################################################################################
################################################################################
################################################################################

## now, after blast, load the contigs that are RdRp positive, and make a subfile of
## those to align with the reference sequences in Geneious

rm(list=ls())

#setwd
homewd <- "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig3/B-RdRp-phylogeny"))

#and import all contigs

#load all positive fecal contigs
fec_contigs_pos <- read.fasta( file = paste0(homewd, "/Fig3/B-RdRp-phylogeny/blast-pos-contigs/all_CoV_pos_contigs_feces.fasta"), as.string = T, forceDNAtolower = F)

#get the names
names_fec <- names(fec_contigs_pos)

#load the names of all the RdRp positive contigs (the hiqual)
RdRP_pos <- read.delim(file = paste0(homewd, "/Fig3/B-RdRP-phylogeny/blast-pos-contigs/20210808_Mada_Bat_CoV_RdRp_blast_feces_nt_results_100len5eval.txt"), header=F)

RdRP_pos <- c(unlist(c(RdRP_pos))) #135
names(RdRP_pos) <- c()

#now subset these names to just get the beginning
RdRP_pos <- sapply(strsplit(RdRP_pos, " "), '[',1)

#now select those sequences from the list (same function as above)
select.seq <- function(seq1, all_pos_contigs){
  
  name_contig = attr(seq1, "name")
  
  #now select that sequence if it matches
  if(length(intersect(name_contig, all_pos_contigs))>0){
    #and return it
    return(seq1)
  }
  
  #otherwise, leave off
  
}

positive_contigs <- lapply(fec_contigs_pos,select.seq, all_pos_contigs=RdRP_pos)

#remove NULL entries
positive_contigs <- positive_contigs[lengths(positive_contigs) > 0] #down to 18 contigs

#and save
names_pos = names(positive_contigs)
write.fasta(positive_contigs, names = names_pos, file.out = paste0(homewd, "/Fig3/B-RdRP-phylogeny/blast-pos-contigs/RdRp_pos_contigs_hiqual_feces.fasta"), as.string = T, open="w")

#now align these in Genious with those used in the database
