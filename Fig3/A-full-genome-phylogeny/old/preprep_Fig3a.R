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



#and call in the AlphaCoVs
#do the same
alpha.dat = read.csv(file = "fullgenomeAlphaCoVs.csv", header = T, stringsAsFactors = F)
head(alpha.dat)
sort(unique(alpha.dat$Host))
unique(alpha.dat$Species[alpha.dat$Host==""]) #4 bat CoVs in here

add.on = subset(alpha.dat, Host=="" & Species=="Miniopterus bat coronavirus 1" | 
                           Host=="" & Species=="Scotophilus bat coronavirus 512" | 
                           Host=="" & Species=="Rhinolophus bat coronavirus HKU2" | 
                           Host=="" & Species=="Miniopterus bat coronavirus HKU8" )

bat.alpha = subset(alpha.dat, Host=="Chiroptera" | Host=="Chaerephon plicatus" | Host=="Cynopterus sphinx" |
                              Host=="Desmodus rotundus" | Host=="Hipposideros" | Host=="Hipposideros cineraceus" | Host=="Hipposideros larvatus" | Host=="Hipposideros pomona" |
                              Host=="Macronycteris vittata" | Host=="Miniopterus fuliginosus" | Host=="Miniopterus pusillus" | Host=="Miniopterus schreibersii" |
                              Host=="Murina cyclotis" | Host=="Myotis brandtii" | Host=="Myotis laniger" | Host=="Myotis lucifugus" | Host=="Myotis muricola" |
                              Host=="Myotis ricketti"| Host=="Myotis sp."| Host=="Nyctalus velutinus"| Host=="Pipistrellus kuhlii"| Host=="Rhinolophus affinis"| Host=="Rhinolophus ferrumequinum" |
                              Host == "Rhinolophus malayanus" | Host=="Rhinolophus sinicus" | Host=="Rhinolophus stheno" | Host=="Rousettus aegyptiacus" | Host=="Scotophilus kuhlii" | Host=="Triaenops afer" | Host=="Tylonycteris robustula")


nrow(bat.alpha) + nrow(add.on) #118 bat full genome alphaCoVs

98+118 #216 full genome bat CoVs

#now get refs
ref.alpha = subset(alpha.dat, Sequence_Type=="RefSeq")
sort(unique(ref.alpha$Host))
#keep all those without bat hosts

ref.alpha = subset(ref.alpha, Host=="" | Host=="Apodemus chevrieri" | Host=="Camelus" | Host=="Mustela putorius" |
                              Host=="Neovison vison" | Host == "Rattus norvegicus" | Host =="Suncus murinus" | Host=="Sus scrofa")

ref.alpha$Species[ref.alpha$Host==""]
alpha.CoV <- rbind(bat.alpha, add.on, ref.alpha)
alpha.CoV <- alpha.CoV[!duplicated(alpha.CoV),] #removes 4 found in the bat data above


#now add in those that are unclassified and have bat hosts
unclass.CoV <- read.csv(file = "fullgenomeUnclassifiedCoVs.csv", header = T, stringsAsFactors = F)
head(unclass.CoV)
unclass.CoV <- dplyr::select(unclass.CoV, names(alpha.CoV))



#118 bat alphaCoVs
#13 ref alphaCoV that are not bat

#98 bat betaCoV
#14 ref betaCoV that are not bat

#44 bat Unclassified CoVs


#and get the text to download from NCBI
all.CoV <- rbind(betaCoV, alpha.CoV, unclass.CoV) #287 genomes + 3 madagascar sequences
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
