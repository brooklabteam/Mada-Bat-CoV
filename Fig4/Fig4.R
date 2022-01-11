rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
#library(rBt)
library(treeio)


#make Bayesian timetree from Nobecovirus relaxed molecular clock model

#first, read in the tree

homewd= "/Users/carabrook/Developer/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig4"))

tree <- read.beast(file = paste0(homewd, "/Fig4/nobecoRecombAVG.tree"))


treedat <- cbind.data.frame(tip_name = tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name
#tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoStrict/AvgNobecoStrictNexus.trees"))

#tree$node.label <- round(tree$posterior,2)
#treedat <- cbind.data.frame(tip_name = tree$tip.label)
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_009021", "NC_048212")
#names(treedat)[names(treedat)=="tip_name"] <- "beast_name"


#and load data of corresponding tree

dat <- read.csv(file = "fig4beast_nobecov_metadata_manual.csv", header = T, stringsAsFactors = F)
dat$collection_date <- as.Date(dat$collection_date)


#test 

mrsd.dat <- max(dat$collection_date)
p1 <- ggtree(tree, mrsd=mrsd.dat)  + theme_tree2()  +geom_nodelab()

tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

#and 
head(dat)

dat$clade <- dat$strain
dat$clade[is.na(dat$clade)] <- "African Eidolon"
dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 
dat$clade[dat$clade == "GX2018"] <- "BtCoV92 / GX2018" 
dat$clade[dat$clade == "BtCoV92"] <- "BtCoV92 / GX2018" 
#dat$clade[dat$accession_num == "KU182962"] <- "BtCoV92 / GX2018" 
#dat$clade[dat$host == "Eidolon_helvum" | dat$host == "Rousettus_madagascariensis"] <- "African Eidolon" 
#dat$clade[dat$accession_num=="MG693170"] <- "HKU9"
#dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

head(dat.plot)
dat.plot$new_label = NA
dat.plot$new_label[!is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[!is.na(dat.plot$strain)], " | ", 
                                             dat.plot$strain[!is.na(dat.plot$strain)], " | ", 
                                             dat.plot$host[!is.na(dat.plot$strain)], " | ",
                                             dat.plot$country[!is.na(dat.plot$strain)], " | ",
                                             dat.plot$collection_year[!is.na(dat.plot$strain)])

dat.plot$new_label[is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[is.na(dat.plot$strain)], " | ", 
                                          dat.plot$host[is.na(dat.plot$strain)], " | ",
                                          dat.plot$country[is.na(dat.plot$strain)], " | ",
                                          dat.plot$collection_year[is.na(dat.plot$strain)])


tree@phylo$tip.label <- dat.plot$new_label

dat.sub <- dplyr::select(dat.plot, new_label, collection_date, country, clade)
head(dat.sub)
dat.sub$clade <- as.factor(dat.sub$clade)

p2 <-ggtree(tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade)) +
      geom_tiplab(size=3) + geom_nodelab(size=2,nudge_x = -15, nudge_y = .7) +
      theme_tree2() +
      theme(legend.position = c(.1,.85),
            plot.margin = unit(c(2,20,2,3), "lines")) +
      coord_cartesian(clip = "off")

#and a second node label that is the date for the P ruf and the original
orig.date <- round(node.sub$nodetime[33],0)

nodePruf <- MRCA(tree, which(tree@phylo$tip.label == "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018"),which(tree@phylo$tip.label == "MG693172  |  Eidolon_helvum  |  Cameroon  |  2013"))
nodeRous <- MRCA(tree, which(tree@phylo$tip.label == "OK067320  |  Rousettus_madagascariensis  |  Madagascar  |  2018"),which(tree@phylo$tip.label == "HM211099  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005"))
#nodeall <- MRCA(tree, which(tree$tip.label == "KP696747  |  Pteropus_rufus  |  Madagascar  |  2011"),which(tree$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016"))
#orig.date <- round(node.sub$nodetime[nodeall],0)
Pruf.date <- round(node.sub$nodetime[nodePruf],0) #date that P ruf branches off from everything else
recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365
Pruf.mean <- round(recent.age-tree@data$height[35],0) #sane as above--date that P. ruf branches off -
#this is taking the P. ruf tip and subtracting the branch length from present to get its position
Pruf.uci <- round(recent.age-tree@data$height_0.95_HPD[35][[1]][1],0)
Pruf.lci <- round(recent.age-tree@data$height_0.95_HPD[35][[1]][2],0)
#should be close to the lengths: tree@data$length_0.95_HPD[35] (does not have it)
Pruf.date <- paste0("~", Pruf.mean, "\n[", Pruf.lci, "-", Pruf.uci, "]")#from FigTree

Rous.mean <- round(recent.age-tree@data$height[33],0)
Rous.date <- round(node.sub$nodetime[nodeRous],0)#they are the same--we're good
Rous.uci <- round(recent.age-tree@data$height_0.95_HPD[33][[1]][1],0)
Rous.lci <- round(recent.age-tree@data$height_0.95_HPD[33][[1]][2],0)
Rous.date <- paste0("~", Rous.mean, "\n[", Rous.lci, "-", Rous.uci, "]")#from FigTree


new.nodel.lab <- rep(NA, nrow(node.sub))
#new.nodel.lab[nodeall] <- paste0("~",orig.date)
new.nodel.lab[nodePruf] <- Pruf.date
new.nodel.lab[nodeRous] <- Rous.date

dat.sub$clade <- as.character(dat.sub$clade)
dat.sub$clade[dat.sub$clade=="African Eidolon"] <- "African~italic(Eidolon)"
dat.sub$clade[dat.sub$clade=="Madagascar Pteropus"] <- "Madagascar~italic(Pteropus)"

dat.sub$novel = "no"
dat.sub$novel[dat.sub$country=="Madagascar"] <- "yes"

colz2 = c('yes' =  "yellow", 'no' = "white")

p3 <-ggtree(tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade), size=3) +
  #geom_tiplab(size=3, nudge_x=5) + 
  geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  geom_nodelab(aes(label=new.nodel.lab), size=4,nudge_x = -55, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  theme_tree2() +
  #geom_treescale(fontsize=3, x=1300,y=22, linesize = .5, width=200,label="years") + 
  scale_color_discrete(labels=c(parse(text="African~italic(Eidolon)"), "BtCoV92 / GX2018", "HKU9", parse(text="Madagascar~italic(Pteropus)"))) +
  theme(legend.position = c(.2,.6), 
        plot.margin = unit(c(.2,20,2,3), "lines")) +
  coord_cartesian(clip = "off", xlim=c(1600, 2150)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(1700, 1800, 1900, 2000))+ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.05) + scale_fill_manual(values=colz2)
  


ggsave(file = paste0(homewd, "/final-figures/Fig4.png"),
       units="mm",  
       width=90, 
       height=60, 
       #limitsize = F,
       scale=3)#, 


#and the version that includes  the recombinant clade

tree <- read.beast(file = paste0(homewd, "/Fig4/nobecoAllAVG.tree"))


treedat <- cbind.data.frame(tip_name = tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name
#tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoStrict/AvgNobecoStrictNexus.trees"))

#tree$node.label <- round(tree$posterior,2)
#treedat <- cbind.data.frame(tip_name = tree$tip.label)
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_009021", "NC_030886", "NC_048212")
#names(treedat)[names(treedat)=="tip_name"] <- "beast_name"


#and load data of corresponding tree

dat <- read.csv(file = "fig4beast_nobecov_metadata_manual.csv", header = T, stringsAsFactors = F)
dat$collection_date <- as.Date(dat$collection_date)


#test 

mrsd.dat <- max(dat$collection_date)
p1 <- ggtree(tree, mrsd=mrsd.dat)  + theme_tree2()  +geom_nodelab()

tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

#and 
head(dat)

dat$clade <- dat$strain
dat$clade[is.na(dat$clade)] <- "African Eidolon"
dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 
dat$clade[dat$clade == "GX2018"] <- "BtCoV92 / GX2018" 
dat$clade[dat$clade == "BtCoV92"] <- "BtCoV92 / GX2018" 
#dat$clade[dat$accession_num == "KU182962"] <- "BtCoV92 / GX2018" 
#dat$clade[dat$host == "Eidolon_helvum" | dat$host == "Rousettus_madagascariensis"] <- "African Eidolon" 
#dat$clade[dat$accession_num=="MG693170"] <- "HKU9"
#dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

head(dat.plot)
dat.plot$new_label = NA
dat.plot$new_label[!is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[!is.na(dat.plot$strain)], " | ", 
                                                     dat.plot$strain[!is.na(dat.plot$strain)], " | ", 
                                                     dat.plot$host[!is.na(dat.plot$strain)], " | ",
                                                     dat.plot$country[!is.na(dat.plot$strain)], " | ",
                                                     dat.plot$collection_year[!is.na(dat.plot$strain)])

dat.plot$new_label[is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[is.na(dat.plot$strain)], " | ", 
                                                    dat.plot$host[is.na(dat.plot$strain)], " | ",
                                                    dat.plot$country[is.na(dat.plot$strain)], " | ",
                                                    dat.plot$collection_year[is.na(dat.plot$strain)])


tree@phylo$tip.label <- dat.plot$new_label

dat.sub <- dplyr::select(dat.plot, new_label, collection_date, country, clade)
head(dat.sub)
dat.sub$clade <- as.factor(dat.sub$clade)

p2 <-ggtree(tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade)) +
  geom_tiplab(size=3) + geom_nodelab(size=2,nudge_x = -15, nudge_y = .7) +
  theme_tree2() +
  theme(legend.position = c(.1,.85),
        plot.margin = unit(c(2,20,2,3), "lines")) +
  coord_cartesian(clip = "off")

#and a second node label that is the date for the P ruf and the original
orig.date <- round(node.sub$nodetime[33],0)

nodePruf <- MRCA(tree, which(tree@phylo$tip.label == "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018"),which(tree@phylo$tip.label == "MG693172  |  Eidolon_helvum  |  Cameroon  |  2013"))
nodeRous <- MRCA(tree, which(tree@phylo$tip.label == "OK067320  |  Rousettus_madagascariensis  |  Madagascar  |  2018"),which(tree@phylo$tip.label == "HM211099  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005"))
#nodeall <- MRCA(tree, which(tree$tip.label == "KP696747  |  Pteropus_rufus  |  Madagascar  |  2011"),which(tree$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016"))
#orig.date <- round(node.sub$nodetime[nodeall],0)
Pruf.date <- round(node.sub$nodetime[nodePruf],0)
recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365
Pruf.mean <- round(recent.age-tree@data$height[which(tree@data$height==sort(tree@data$height)[length(tree@data$height)])],0)
Pruf.uci <- round(recent.age-tree@data$height_range[which(tree@data$height==sort(tree@data$height)[length(tree@data$height)])][[1]][1],0)
Pruf.lci <- round(recent.age-tree@data$height_range[which(tree@data$height==sort(tree@data$height)[length(tree@data$height)])][[1]][2],0)
Pruf.date <- paste0("~", Pruf.mean, "\n[", Pruf.lci, "-", Pruf.uci, "]")#from FigTree


Rous.mean <- round(recent.age-tree@data$height[which(tree@data$height==sort(tree@data$height)[length(tree@data$height)-1])],0)
Rous.date <- round(node.sub$nodetime[nodeRous],0)
Rous.uci <- round(recent.age-tree@data$height_range[which(tree@data$height==sort(tree@data$height)[length(tree@data$height)-1])][[1]][1],0)
Rous.lci <- round(recent.age-tree@data$height_range[which(tree@data$height==sort(tree@data$height)[length(tree@data$height)-1])][[1]][2],0)
Rous.date <- paste0("~", Rous.mean, "\n[", Rous.lci, "-", Rous.uci, "]")#from FigTree


new.nodel.lab <- rep(NA, nrow(node.sub))
#new.nodel.lab[nodeall] <- paste0("~",orig.date)
new.nodel.lab[nodePruf] <- Pruf.date
new.nodel.lab[nodeRous] <- Rous.date

dat.sub$clade <- as.character(dat.sub$clade)
dat.sub$clade[dat.sub$clade=="African Eidolon"] <- "African~italic(Eidolon)"
dat.sub$clade[dat.sub$clade=="Madagascar Pteropus"] <- "Madagascar~italic(Pteropus)"
dat.sub$novel = "no"
dat.sub$novel[dat.sub$country=="Madagascar"] <- "yes"

colz2 = c('yes' =  "yellow", 'no' = "white")



p3 <-ggtree(tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade), size=3) +
  #geom_tiplab(size=3, nudge_x=5) + 
  geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  geom_nodelab(aes(label=new.nodel.lab), size=4,nudge_x = -50, nudge_y = -1,  color="firebrick", fontface=2, geom="label", fill="white") +
  theme_tree2() +
  #geom_treescale(fontsize=3, x=1300,y=22, linesize = .5, width=200,label="years") + 
  scale_color_discrete(labels=c(parse(text="African~italic(Eidolon)"), "BtCoV92 / GX2018", "GCCDC1", "HKU9", parse(text="Madagascar~italic(Pteropus)"))) +
  theme(legend.position = c(.18,.7), 
        plot.margin = unit(c(.2,20,2,3), "lines")) +
  coord_cartesian(clip = "off", xlim=c(1600, 2150)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(1700, 1800, 1900, 2000)) +ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.05) + scale_fill_manual(values=colz2)



ggsave(file = paste0(homewd, "/final-figures/FigS3.png"),
       units="mm",  
       width=90, 
       height=60, 
       #limitsize = F,
       scale=3)#, 

