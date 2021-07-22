#rename your fasta files to something shorter before adding the filename to all the headers

rm(list=ls())

library(seqinr)

#now make a list of the files in this directory (each is a fasta except the R file)
all_names <- c(unlist(list.files()))

#all the names except this R file
all_names <- all_names[all_names!="rename-fastas.R"]

#now make new names
new_names <- paste0(sapply(strsplit(all_names, split="_fecRa"), '[',1),".fasta")
new_names <- as.list(new_names)
all_names <- as.list(all_names)
#now load and rename

rename.func <- function(old, new){
  input = read.fasta(old, as.string=T, forceDNAtolower = F)
  names_input = names(input)
  names_input = strsplit(names_input, split = "_cov")
  names_input = c(unlist(sapply(names_input, '[', 1)))
  new_slim = sapply(strsplit(old, split="_fecRa"), '[',1)
  
  names_output = paste(new_slim, names_input, sep="_")
  
  new_name = paste0(new_slim, ".fasta")
  
  write.fasta(input, names=names_output, file.out = new_name)
  
}



