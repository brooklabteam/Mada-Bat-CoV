rm(list=ls())

library(plyr)
library(dplyr)


#set wd
homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig3/A-full-genome-phylogeny/"))

#load the dataset and query
dat <- read.csv(file = "fullgenomeBetaCoVs.csv", header = T, stringsAsFactors = F)

head(dat)
#look at the unique hosts
sort(unique(dat$Host))
sort(unique(dat$Species[dat$Host==""])) #2 bat viruses
add.on.beta = subset(dat, Species=="Pipistrellus bat coronavirus HKU5" & Host=="" | Species=="Tylonycteris bat coronavirus HKU4" & Host=="")

#select all the bats
bat.dat = subset(dat, Host=="Chiroptera" | Host =="Aselliscus stoliczkanus" | Host =="Chaerephon plicatus"|
                      Host=="Chaerephon plicatus" | Host == "Eonycteris spelaea" | Host =="Hipposideros pratti" | Host =="Hypsugo pulveratus"|
                      Host== "Hypsugo savii" | Host == "Laephotis capensis" | Host =="Pipistrellus abramus" | Host =="Pipistrellus kuhlii" |
                      Host == "Rhinolophus" | Host == "Rhinolophus affinis" | Host == "Rhinolophus blasii" | Host == "Rhinolophus cornutus" |
                      Host == "Rhinolophus ferrumequinum" | Host == "Rhinolophus hipposideros" | Host == "Rhinolophus macrotis" |
                      Host == "Rhinolophus malayanus" | Host == "Rhinolophus pusillus"  | Host == "Rhinolophus sinicus" | Host == "Rhinolophus stheno" |
                      Host == "Rousettus leschenaultii" | Host == "Rousettus sp."  | Host == "Tylonycteris pachypus" | Host == "Vespertilio sinensis")

sort(unique(bat.dat$Species))

nrow(bat.dat) + nrow(add.on.beta) #98 bat betaCoV full genomes

#then select the ref seq for the hosts that are not bats
ref.dat = subset(dat, Sequence_Type=="RefSeq")
sort(unique(ref.dat$Host))
subset(ref.dat, Host=="")#check those with no host... none are batCoVs

#now take just those with no bats...
ref.sub = subset(ref.dat, Host!="Chiroptera" & Host!="Rousettus leschenaultii" & Host !="Rhinolophus blasii")

#now put these together to draw from GenBank
betaCoV <- rbind(bat.dat, add.on.beta, ref.sub)

betaCoV <- betaCoV[!duplicated(betaCoV),] #removes one of the ref found in the bat

#98 bat betaCoV
#14 ref betaCoV that are not bat


#now add in those that are unclassified betaCoVs and have bat hosts

unclass.CoV<- read.csv(file = "fullgenomeUnclassifiedCoVs_manual.csv", header = T, stringsAsFactors = F)
head(unclass.CoV)
names(unclass.CoV)[names(unclass.CoV)=="Accession_Number"] <- "Accession"
#and take only those that are betaCoVs:
unique(unclass.CoV$sub_group)
unclass.CoV = subset(unclass.CoV, sub_group!="Alphacoronavirus")
unclass.CoV <- dplyr::select(unclass.CoV, names(betaCoV))

unclass.CoV <- unclass.CoV[!duplicated(unclass.CoV),] #removes 6 duplicates

#98 bat betaCoV
#14 ref betaCoV that are not bat

#24 bat Unclassified betaCoVs


#and get the text to download from NCBI
all.CoV <- rbind(betaCoV, unclass.CoV) #136 genomes + 3 madagascar sequences
accession_num <- paste(c(all.CoV$Accession), collapse = ", ")

#now put this into your webbrowser to download
text.for.NCBI <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=",accession_num)

#once downloaded, send to MAFFT for alignment
rm(list=ls())
#then, after alignment is ready, prepare the names for RAxML (no space, semicolon, colon, parentheses, dash, slash, comma, quote allowed in name (should just all be underscore)
library(seqinr)
#library(msa)
alignment1 <- read.alignment(file = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/Fig3/A-full-genome-phylogeny/allCoVsalign825.fasta", format="fasta", forceToLower = F) #alignment_fullgenomeCoVs_8_7_rename.fasta", format="fasta", forceToLower = F)

                             
tmp <- as.list(alignment1$nam)

change.spacing <- function(df){
  df_new <- sapply(strsplit(df,"-"), function(x) x[[1]])
  return(df_new)
}

names_new = c(unlist(lapply(tmp, change.spacing)))
#alignment1$seq[[1]]

#new_names <- sub("__", "_", new_names) 
class(alignment1$seq)
write.fasta(sequences = as.list(alignment1$seq), names = names_new, file.out =  "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/Fig3/A-full-genome-phylogeny/RAxML_alignment_fullgenomeCoVs_8_25.fasta", as.string = T, open="w")

#now send aboce to modeltest and eventually RAxML
