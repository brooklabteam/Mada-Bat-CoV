rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(seqinr)
library(ggmsa)


#load the msa of just the nobecovirus genomes



homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"
setwd(paste0(homewd, "Fig2B")) 

msa = tidy_msa(msa="AllNobecoSubAlign.fasta", start = 28449, end = 31558)
head(msa)
msa$clade <- NA
msa$clade[msa$name=="OK067319_Pteropus_rufus_AMB130"] <- "Madagascar Pteropus"
msa$clade[msa$name=="MK211379_1_GX2018_Rhinolophus_affinis_China"] <- "BtCoV92/GX2018"
msa$clade[msa$name=="MK492263_1_BtCoV92_Cynopteris_brachyotis_Singapore"] <- "BtCoV92/GX2018"
msa$clade[msa$name=="NC_030886_GCCDC1_Rousettus_leschenaulti_China"] <- "GCCDC1"
msa$clade[msa$name=="MT350598_1_GCCDC1_Eonycteris_spelaea_Singapore"] <- "GCCDC1"
msa$clade[msa$name=="OK067321_Rousettus_madagascariensis_MIZ240"] <- "African Eidolon"
msa$clade[msa$name=="OK067320_Rousettus_madagascariensis_MIZ178"] <- "African Eidolon"
msa$clade[msa$name=="NC_048212_1_Eidolon_helvum_Cameroon"] <- "African Eidolon"
msa$clade[msa$name=="NC_009021_1_HKU9_Rousettus_leschenaulti_China"] <- "HKU9"
msa$clade[msa$name=="MG762674_1_HKU9_Rousettus_sp_China"] <- "HKU9"

names(msa)[names(msa)=="character"] <- "nucleotide"
head(msa)

#get new names
dat <- read.csv(file = paste0(homewd, "Fig2B/fig2B_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)[names(dat)=="tip_label"] <- "name"

dat$name_new[!is.na(dat$strain)] <-paste(dat$accession_num[!is.na(dat$strain)], " | ",
                                            dat$strain[!is.na(dat$strain)], " | ",
                                            dat$host[!is.na(dat$strain)], " | ", 
                                            dat$country[!is.na(dat$strain)], " | ",
                                            dat$collection_year[!is.na(dat$strain)])

dat$name_new[is.na(dat$strain)] <-paste(dat$accession_num[is.na(dat$strain)], " | ",
                                           dat$host[is.na(dat$strain)], " | ", 
                                           dat$country[is.na(dat$strain)], " | ",
                                           dat$collection_year[is.na(dat$strain)])


dat <- dplyr::select(dat, name, name_new)
msa.dat <- merge(msa, dat, by="name", all.x=T, sort=F)

unique(msa.dat$name_new)
#and the border
border.dat <- data.frame(name=unique(msa.dat$name_new), position=rep(c(28455,30250),each=10))
border.dat$name <- factor(border.dat$name, levels=c( "OK067321  |  Rousettus_madagascariensis  |  Madagascar  |  2018", "OK067320  |  Rousettus_madagascariensis  |  Madagascar  |  2018", "NC_048212  |  Eidolon_helvum  |  Cameroon  |  2013", "NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014", "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016" , "NC_009021  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005", "MG762674  |  HKU9  |  Rousettus_sp  |  China  |  2009", "MK492263  |  BtCoV92  |  Cynopteris_brachyotis  |  Singapore  |  2015", "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016", "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018"))

msa.dat$name <- msa.dat$name_new
#add genome details
genome.df <- data.frame(name= "genome annotation",
                        start = c(28449, 28846, 29243), stop=c(28846, 29243,30263), 
                        gene = c("N","p10", "NS7"))
genome.df$labx <- genome.df$stop - (genome.df$stop-genome.df$start)/2
genome.df

highlight.df <- data.frame(name= rep(c("NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014", "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016"), each=2),
                           position = rep(c(28846, 29243), 2))

highlight.df <- data.frame(name= c("NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014", "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016"),
                           start = rep(c(28846),2), stop= rep(29243, 2))


#msa.all <- rbind(msa, genome.df)

#msa.all$name <- factor(msa.all$name, levels=c("genome annotation", "OK067321_Rousettus_madagascariensis_MIZ240", "OK067320_Rousettus_madagascariensis_MIZ178", "NC_048212_1_Eidolon_helvum_Cameroon", "NC_030886_GCCDC1_Rousettus_leschenaulti_China", "MT350598_1_GCCDC1_Eonycteris_spelaea_Singapore", "NC_009021_1_HKU9_Rousettus_leschenaulti_China", "MG762674_1_HKU9_Rousettus_sp_China", "MK492263_1_BtCoV92_Cynopteris_brachyotis_Singapore", "MK211379_1_GX2018_Rhinolophus_affinis_China", "OK067319_Pteropus_rufus_AMB130"))
msa.dat$name <- factor(msa.dat$name, levels=c( "OK067321  |  Rousettus_madagascariensis  |  Madagascar  |  2018", "OK067320  |  Rousettus_madagascariensis  |  Madagascar  |  2018", "NC_048212  |  Eidolon_helvum  |  Cameroon  |  2013", "NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014", "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016" , "NC_009021  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005", "MG762674  |  HKU9  |  Rousettus_sp  |  China  |  2009", "MK492263  |  BtCoV92  |  Cynopteris_brachyotis  |  Singapore  |  2015", "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016", "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018"))

#msa.dat$clade[msa.dat$clade=="African Eidolon"] <- 
colz = c("A" = "seagreen3", "T" = "cornflowerblue", "G"="tomato", "C" = "plum3", "-" = "white")# "N"= "mediumpurple1", "p10"="red", "NS7" = "magenta")

msa.plot <- ggplot(data=subset(msa.dat, position >= 28449 & position <=30263)) + 
                geom_tile(aes(x=position, y=name, fill=nucleotide, color=nucleotide),width=0.9, height=0.9) +
                geom_rect(data=highlight.df, aes(xmin=start, xmax=stop, ymin=3.5, ymax=4.5), color="orange", fill="gray", size=3, alpha=.2) +
                geom_rect(data=highlight.df, aes(xmin=start, xmax=stop, ymin=4.5, ymax=5.5), color="orange", fill="gray", size=3, alpha=.2) +
                scale_fill_manual(values=colz) + theme_bw() + 
                scale_color_manual(values=colz) +
                theme(axis.title.y = element_blank(), 
                      panel.background = element_rect(fill="black"),
                      #panel.border = element_line(color="black", size=1),
                      panel.grid = element_blank(),
                      axis.text = element_text(size=14),
                      axis.ticks = element_blank(),
                      axis.title.x = element_text(size=16),
                      legend.text = element_text(size=10),
                      plot.margin = unit(c(.1, 2,.1, .1), "cm"),
                      legend.position = "bottom", legend.direction = "horizontal",legend.title=element_blank()) +
                xlab("genome position")+ coord_cartesian(expand = F, xlim=c(28449, 30263), clip = "off")  +
        annotate("text", x=30300, y=10, label="atop(Madagascar,italic(Pteropus))", angle=270, size=5, parse=T) +
        annotate("text", x=30300, y=8.5, label="BtCoV92 / GX2018", angle=270, size=5) +
        annotate("text", x=30300, y=6.5, label="HKU9", angle=270, size=5) +
        annotate("text", x=30300, y=4.5, label="GCCDC1", angle=270, size=5) +
        annotate("text", x=30300, y=2, label="African~italic(Eidolon)", angle=270, size=5, parse=T) 

genome.plot <- ggplot(data=genome.df) + 
        geom_rect(aes(xmin=start, xmax=stop, ymin=0, ymax=1), fill="white", color="black", size=2) +
        theme_bw() + 
        theme(panel.grid = element_blank(), axis.title = element_blank(),
              axis.ticks = element_blank(), axis.text=element_blank(),
              plot.margin = unit(c(.1, 2,0, 15.5), "cm")) +
        coord_cartesian(expand = F, xlim=c(28449, 30263)) +
        geom_label(aes(x=labx, y=.5, label=gene), label.size = NA, size=5)
        

msa.all <- cowplot::plot_grid(genome.plot, msa.plot, nrow=2, ncol=1, rel_heights = c(.06,1))


ggsave(file = paste0(homewd, "/final-figures/Fig2B.png"),
       plot = msa.all,
       units="mm",  
       width=180, 
       height=100, 
       #limitsize = F,
       scale=3)#, 


