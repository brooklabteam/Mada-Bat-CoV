rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)
library(rentrez) #use this program to query records from NCBI

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"

setwd(paste0(homewd, "/Fig3"))

#load the fig3a tree
treeA <-  read.tree(file = paste0(homewd, "Fig3/A-full-genome-phylogeny/Fig3A-raxml-output/Fig3A.raxml.supportTBE"))

#get the accession numbers
dat <- treeA$tip.label

accession_nums <- paste(lapply( strsplit(dat, "_"), function(x) x[1]), lapply( strsplit(dat, "_"), function(x) x[2]), sep="_")

#and those that are new:
accession_nums = subset(accession_nums, accession_nums!="Bat_CoV" & accession_nums!="BatCoV_R")
accession_nums[accession_nums=="MN996532_2"] <- "MN996532"
accession_nums[accession_nums=="KC869678_4"] <- "KC869678"

accession_nums <- sub("_1", "", accession_nums)


#now write a function to query all accession numbers and collect
#relevant meta data
query.meta <- function(acc_num1){
   out <- entrez_summary(db="nuccore",
                 id=acc_num1,
                 version = "2.0")
   
   
   out2 <- strsplit(out$subname, "|", fixed = T)
   out3 <- strsplit(out$subtype, "|", fixed = T)
   
   out.df <- rbind.data.frame(unlist(out2))
   names(out.df) <- rbind(unlist(out3))
   
   #select only those you want 
   new.df <- cbind.data.frame(strain="NA", host="NA", country="NA", collection_date="NA")
   new.df$strain <- if(length(out.df$strain)>0){out.df$strain}else{"NA"}
   new.df$host <- if(length(out.df$host)>0){out.df$host}else{"NA"}
   new.df$country <- if(length(out.df$country)>0){out.df$country}else{"NA"}
   new.df$collection_date <- if(length(out.df$collection_date)>0){out.df$collection_date}else{"NA"}
   
   
   
   #and add
   new.df$accession_num = acc_num1
   new.df$title = if(length(out$title)>0){out$title}else{"NA"}
   new.df$sub_group <-if(length(out$organism)>0){out$organism}else{"NA"}
   
   #remove
   new.df$title <- sub(", complete genome", "", new.df$title)
   
   #and add in sub-genus manually
   return(new.df)
   
}


accession_num_list <- as.list(accession_nums)

#and apply function
accession_meta_list <- lapply(accession_num_list, query.meta)

accession.df <- data.table::rbindlist(accession_meta_list)
head(accession.df)
#and link to the tip names on the tree
dat <- data.frame(tip_label= treeA$tip.label)

dat$accession_num <- paste(lapply( strsplit(dat$tip_label, "_"), function(x) x[1]), lapply( strsplit(dat$tip_label, "_"), function(x) x[2]), sep="_")
head(dat)

#and those that are new:
dat = subset(dat, accession_num!="Bat_CoV" & accession_num!="BatCoV_R")
dat$accession_num[dat$accession_num=="MN996532_2"] <- "MN996532"
dat$accession_num[dat$accession_num=="MN996532_2"] <- "MN996532"
dat$accession_num[dat$accession_num=="KC869678_4"] <- "KC869678"

dat$accession_num <- sub("_1", "", dat$accession_num)
accession.df <- merge(accession.df, dat, by="accession_num", all.x=T)
head(accession.df)
write.csv(accession.df, file = paste0(homewd,"/Fig3/A-full-genome-phylogeny/fig3a_metadata.csv"), row.names=F)

#now curate manually and re-import.

############################
rm(list=ls())
###panel B

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"

setwd(paste0(homewd, "/Fig3"))

treeB <-  read.tree(file = paste0(homewd, "Fig3/B-RdRP-phylogeny/Fig3B-raxml-output/Fig3B.raxml.supportTBE"))


#take a quick look in base R
plot(treeB)

#get the accession numbers
dat <- treeB$tip.label

accession_nums <- paste(lapply( strsplit(dat, "_"), function(x) x[3]), lapply( strsplit(dat, "_"), function(x) x[4]), sep="_")


#and those that are new:
accession_nums = subset(accession_nums, 
                        accession_nums!="RR034B_220" & 
                        accession_nums!="RR034B_288" & 
                        accession_nums!="RR034B_232" &
                           accession_nums!="RR034B_244" &
                           accession_nums!="RR034B79" & 
                           accession_nums!="RR034B_021" & 
                           accession_nums!="RR034B_010")

accession_nums <- sub("_1", "", accession_nums)

#edit a few
accession_nums[accession_nums=="010800"] <- "NC_010800"
accession_nums[accession_nums=="009020"] <- "NC_009020"
accession_nums[accession_nums=="039207"] <- "NC_039207"


#now write a function to query all accession numbers and collect
#relevant meta data
query.meta <- function(acc_num1){
   out <- entrez_summary(db="nuccore",
                         id=acc_num1,
                         version = "2.0")
   
   
   out2 <- strsplit(out$subname, "|", fixed = T)
   out3 <- strsplit(out$subtype, "|", fixed = T)
   
   out.df <- rbind.data.frame(unlist(out2))
   names(out.df) <- rbind(unlist(out3))
   
   #select only those you want 
   new.df <- cbind.data.frame(strain="NA", host="NA", country="NA", collection_date="NA")
   new.df$strain <- if(length(out.df$strain)>0){out.df$strain}else{"NA"}
   new.df$host <- if(length(out.df$host)>0){out.df$host}else{"NA"}
   new.df$country <- if(length(out.df$country)>0){out.df$country}else{"NA"}
   new.df$collection_date <- if(length(out.df$collection_date)>0){out.df$collection_date}else{"NA"}
   
   
   
   #and add
   new.df$accession_num = acc_num1
   new.df$title = if(length(out$title)>0){out$title}else{"NA"}
   new.df$sub_group <-if(length(out$organism)>0){out$organism}else{"NA"}
   
   #remove
   new.df$title <- sub(", complete genome", "", new.df$title)
   
   #and add in sub-genus manually
   return(new.df)
   
}

accession_num_list <- as.list(accession_nums)

#and apply function
accession_meta_list <- lapply(accession_num_list, query.meta)

accession.df <- data.table::rbindlist(accession_meta_list)
head(accession.df)
#and link to the tip names on the tree
dat <- data.frame(tip_label= treeB$tip.label)

dat$accession_num <- paste(lapply( strsplit(dat$tip_label, "_"), function(x) x[3]), lapply( strsplit(dat$tip_label, "_"), function(x) x[4]), sep="_")
head(dat)

#and those that are new:
dat = subset(dat, accession_nums!="RR034B_220" & 
                accession_nums!="RR034B_288" & 
                accession_nums!="RR034B_232" &
                accession_nums!="RR034B_244" &
                accession_nums!="RR034B79" & 
                accession_nums!="RR034B_021" & 
                accession_nums!="RR034B_010")
dat$accession_num[dat$accession_num=="010800"] <- "NC_010800"
dat$accession_num[dat$accession_num=="009020"] <- "NC_009020"
dat$accession_num[dat$accession_num=="039207"] <- "NC_039207"

dat$accession_num <- sub("_1", "", dat$accession_num)
accession.df <- merge(accession.df, dat, by="accession_num", all.x=T)
head(accession.df)

write.csv(accession.df, file = paste0(homewd,"/Fig3/B-RdRP-phylogeny/fig3b_metadata.csv"), row.names=F)





