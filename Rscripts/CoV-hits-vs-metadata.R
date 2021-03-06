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
all.unique <- ddply(all.hits, .(sampleID, tissue), summarise)

#now load our original dataset and compare those that match as positive.

#load Gwen's metadata for positives from IDseq
fec.dat <- read.csv(file = paste0(homewd,"/metadata/bat_fecal_metadata_April2021.csv"), header = T, stringsAsFactors = F)
ur.dat <- read.csv(file = paste0(homewd,"/metadata/bat_urine_WB_metadata_April2021.csv"), header = T, stringsAsFactors = F)

#and compare the positives
IDseq.pos.fec <- data.frame(sampleID =sapply(strsplit(fec.dat$sample_name[fec.dat$CoV=="Y"], "_fec"), '[',1), tissue="feces")
IDseq.pos.ur <- data.frame(sampleID =sapply(strsplit(ur.dat$sample_name[ur.dat$CoV=="Y"], "_ua"), '[',1), tissue="urine")

IDseq.unique <- rbind(IDseq.pos.fec, IDseq.pos.ur)

#now compare

#those in my pipeline but not on IDseq:
setdiff(all.unique$sampleID[all.unique$tissue=="feces"], IDseq.unique$sampleID[IDseq.unique$tissue=="feces"])
#"RR034B_019". actually does come up on IDseq but as a trace hit. 
#contig is good (216 bp in length), so we will keep it

#and those in IDseq but not my pipeline:
setdiff(IDseq.unique$sampleID[IDseq.unique$tissue=="feces"], all.unique$sampleID[all.unique$tissue=="feces"])

#lots. should manually double-check



#"RR034B_004" - keep, was missed, 18 reads, 2 contigs
#"RR034B_016"-
#"RR034B_057"
#"RR034B_058" 
#"RR034B_060"
#"RR034B_061" 
#"RR034B_062"
#"RR034B_064"
#"RR034B_070"
#"RR034B_074"
#"RR034B_099"
#"RR034B_101"
#"RR034B_111"
#"RR034B_112"
#"RR034B_113" 
#"RR034B_121"
#"RR034B_122"
#"RR034B_126"
#"RR034B_129"
#"RR034B_143"
#"RR034B_145"
#"RR034B_148"
#"RR034B_155"
#"RR034B_156"
#"RR034B_158"
#"RR034B_160"
#"RR034B_162"
#"RR034B_164"
#"RR034B_165"
#"RR034B_167"
#"RR034B_168"
#"RR034B_174"
#"RR034B_175" - 
#"RR034B_176" - 14 reads, maybe 1 contig           revisit
#"RR034B_186" - 0 reads, 0 contigs, discard
#"RR034B_188" - 2 reads, 0 contigs, discard
#"RR034B_218" - 39 reads, maybe 2 contigs               revisit
#"RR034B_220" - 185 reads, 17 contigs, include
#"RR034B_228" - 25 reads, no contig,                    revisit
#"RR034B_239" - 1 contig, 11K reads, great hits, include
#"RR034B_242" - 4 contigs, 96 reads, include
#"RR034B_244" - 26 contigs, 292 reads, must include
#"RR034B_250" - 4 reads, no contigs
#"RR034B_256" - 4 reads, no contigs
#"RR034B_260" - 2 reads only, no contigs - discard
#"RR034B_262" - hits to DeltaCoV, 2 reads, no contigs - discard

#and urine - in mine but not Gwen's
setdiff(all.unique$sampleID[all.unique$tissue=="urine"], IDseq.unique$sampleID[IDseq.unique$tissue=="urine"])
#"RR034B_342" -  4 reads, maybe 1 contig. low cov. discard
#"RR034B_401" - 2 reads, 0 contigs. low cov. discard
# "RR034B_404" - 4 reads, 1 contig. low cov. discard
# "RR034B_421" - 4 reads, 0 contig. low cov. discard
# "RR034B_449" - 10 reads, 1 contig. low cov. discard
# "RR034B_506" - 28 reads, 2 contigs, should be real. keep

#and in Gwen's but not mine
setdiff(IDseq.unique$sampleID[IDseq.unique$tissue=="urine"], all.unique$sampleID[all.unique$tissue=="urine"])
# "RR034B_351"
# "RR034B_358"
# "RR034B_362"
# "RR034B_374"
# "RR034B_377"
# "RR034B_395"
# "RR034B_402"
# "RR034B_465"
# "RR034B_474"
# "RR034B_488"
# "RR034B_495"


#upon manual inspection...