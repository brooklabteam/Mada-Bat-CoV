rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
library(rBt)


#make Bayesian timetree from Nobecovirus relaxed molecular clock model

#first, read in the tree

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig4"))

tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoRelaxed/AvgNobecoRelaxedNexus"))
#tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoStrict/AvgNobecoStrictNexus.trees"))

tree$node.label <- round(tree$posterior,2)
treedat <- cbind.data.frame(tip_name = tree$tip.label)
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "-"), function(x) x[[1]])
names(treedat)[names(treedat)=="tip_name"] <- "beast_name"


#and load data of corresponding tree

dat <- read.csv(file = "fig4_Bayesian_tree_RdRp_metadata.csv", header = T, stringsAsFactors = F)
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
dat$clade[dat$clade == "GX2018"] <- "BtCoV92 / GX2018" 
dat$clade[dat$clade == "BtCoV92"] <- "BtCoV92 / GX2018" 
dat$clade[dat$accession_num == "KU182962"] <- "BtCoV92 / GX2018" 
dat$clade[dat$host == "Eidolon_helvum" | dat$host == "Rousettus_madagascariensis"] <- "African Eidolon" 
dat$clade[dat$accession_num=="MG693170"] <- "HKU9"
dat$clade[dat$host == "Pteropus_rufus" ] <- "Madagascar Pteropus" 

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)



dat.plot$new_label[!is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[!is.na(dat.plot$strain)], " | ", 
                                             dat.plot$strain[!is.na(dat.plot$strain)], " | ", 
                                             dat.plot$host[!is.na(dat.plot$strain)], " | ",
                                             dat.plot$country[!is.na(dat.plot$strain)], " | ",
                                             dat.plot$collection_year[!is.na(dat.plot$strain)])

dat.plot$new_label[is.na(dat.plot$strain)] <- paste(dat.plot$accession_num[is.na(dat.plot$strain)], " | ", 
                                          dat.plot$host[is.na(dat.plot$strain)], " | ",
                                          dat.plot$country[is.na(dat.plot$strain)], " | ",
                                          dat.plot$collection_year[is.na(dat.plot$strain)])


tree$tip.label <- dat.plot$new_label

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

nodePruf <- MRCA(tree, which(tree$tip.label == "KP696747  |  Pteropus_rufus  |  Madagascar  |  2011"),which(tree$tip.label == "OK020087  |  Pteropus_rufus  |  Madagascar  |  2018"))
nodeall <- MRCA(tree, which(tree$tip.label == "KP696747  |  Pteropus_rufus  |  Madagascar  |  2011"),which(tree$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016"))
orig.date <- round(node.sub$nodetime[nodeall],0)
#Pruf.date <- round(node.sub$nodetime[nodePruf],0)
recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365
Pruf.mean <- round(recent.age-164.7363,0)
Pruf.uci <- round(recent.age-19.2614,0)
Pruf.lci <- round(recent.age-513.3403,0)
Pruf.date <- paste0("~", Pruf.mean, "\n[", Pruf.lci, "-", Pruf.uci, "]")#from FigTree
new.nodel.lab <- rep(NA, nrow(node.sub))
#new.nodel.lab[nodeall] <- paste0("~",orig.date)
new.nodel.lab[nodePruf] <- Pruf.date

dat.sub$clade <- as.character(dat.sub$clade)
dat.sub$clade[dat.sub$clade=="African Eidolon"] <- "African~italic(Eidolon)"
dat.sub$clade[dat.sub$clade=="Madagascar Pteropus"] <- "Madagascar~italic(Pteropus)"

p3 <-ggtree(tree, mrsd=mrsd.dat) %<+% dat.sub + geom_tippoint(aes(color=clade), size=3) +
  geom_tiplab(size=3, nudge_x=5) + geom_nodelab(size=2.5,nudge_x = -21, nudge_y = .7) +
  geom_nodelab(aes(label=new.nodel.lab), size=3,nudge_x = -66, nudge_y = -1.1,  color="firebrick", fontface=2, geom="label", fill="white") +
  #theme_tree2() +
  geom_treescale(fontsize=3, x=1300,y=22, linesize = .5, width=200,label="years") + 
  scale_color_discrete(labels=c(parse(text="African~italic(Eidolon)"), "BtCoV92 / GX2018", "GCCDC1", "HKU9", parse(text="Madagascar~italic(Pteropus)"))) +
  theme(legend.position = c(.2,.85), 
        plot.margin = unit(c(2,20,2,3), "lines")) +
  coord_cartesian(clip = "off")



ggsave(file = paste0(homewd, "/final-figures/Fig4.png"),
       units="mm",  
       width=90, 
       height=60, 
       #limitsize = F,
       scale=3)#, 

