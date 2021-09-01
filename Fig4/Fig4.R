rm(list = ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

#load trees and make into multipanel amino acid tree plot

homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"
setwd(paste0(homewd, "Fig4")) 

#load trees
Stree  <- read.tree(file = paste0(homewd,"Fig4/4-raxml-output/S/Sgene.raxml.supportFBP"))
Mtree  <- read.tree(file = paste0(homewd,"Fig4/4-raxml-output/M/Mgene.raxml.supportFBP"))
Etree  <- read.tree(file = paste0(homewd,"Fig4/4-raxml-output/E/Egene.raxml.supportFBP"))
Ntree  <- read.tree(file = paste0(homewd,"Fig4/4-raxml-output/N/Ngene.raxml.supportFBP"))

Stree.root<- root(Stree, which(Stree$tip.label == "GammaCoV_NC_010800_1_Turkey_S"))
Mtree.root<- root(Mtree, which(Mtree$tip.label == "GammaCoV_NC_010800_1_Turkey_M"))
Etree.root<- root(Etree, which(Etree$tip.label == "GammaCoV_NC_010800_1_Turkey_E"))
Ntree.root<- root(Ntree, which(Ntree$tip.label == "GammaCoV_NC_010800_1_Turkey_N"))


#load metadata
meta.dat <- read.csv(file = paste0(homewd,"Fig4/amino_acid_metadata_manual.csv"), header = T, stringsAsFactors = F)
head(meta.dat)

meta.dat$new_label <- NA
meta.dat$new_label[meta.dat$strain!=""] <- paste0(meta.dat$accession_number[meta.dat$strain!=""], " | ", 
                                                  meta.dat$strain[meta.dat$strain!=""], " | ",
                                                  meta.dat$host[meta.dat$strain!=""], " | ",
                                                  meta.dat$country[meta.dat$strain!=""], " | ",
                                                  meta.dat$year[meta.dat$strain!=""])

meta.dat$new_label[meta.dat$strain==""] <- paste0(meta.dat$accession_number[meta.dat$strain==""], " | ", 
                                                  meta.dat$host[meta.dat$strain==""], " | ",
                                                  meta.dat$country[meta.dat$strain==""], " | ",
                                                  meta.dat$year[meta.dat$strain==""])
                                                  
meta.dat$novel = 0
meta.dat$novel[meta.dat$country=="Madagascar"] <- 1
meta.dat$novel <- as.factor(meta.dat$novel)
meta.dat$bat_host[meta.dat$bat_host==1] <- "bat host"
meta.dat$bat_host[meta.dat$bat_host==0] <- "non-bat host"
meta.dat$bat_host <- as.factor(meta.dat$bat_host)
metaS <- meta.dat
metaE <- meta.dat
metaN <- meta.dat
metaM <- meta.dat
metaS$tip_label <- paste0(metaS$tip_label, "_S")
metaE$tip_label <- paste0(metaE$tip_label, "_E")
metaN$tip_label <- paste0(metaN$tip_label, "_N")
metaM$tip_label <- paste0(metaM$tip_label, "_M")

Sdat <- data.frame(tip_label=Stree.root$tip.label)
Sdat <- merge(Sdat, metaS, by="tip_label", all.x = T, sort = F)
Sdat$old_label <- Sdat$tip_label
Sdat$tip_label <- Sdat$new_label
Stree.root$tip.label <- Sdat$tip_label

Edat <- data.frame(tip_label=Etree.root$tip.label)
Edat <- merge(Edat, metaE, by="tip_label", all.x = T, sort = F)
Edat$old_label <- Edat$tip_label
Edat$tip_label <- Edat$new_label
Etree.root$tip.label <- Edat$tip_label

Mdat <- data.frame(tip_label=Mtree.root$tip.label)
Mdat <- merge(Mdat, metaM, by="tip_label", all.x = T, sort = F)
Mdat$old_label <- Mdat$tip_label
Mdat$tip_label <- Mdat$new_label
Mtree.root$tip.label <- Mdat$tip_label


Ndat <- data.frame(tip_label=Ntree.root$tip.label)
Ndat <- merge(Ndat, metaN, by="tip_label", all.x = T, sort = F)
Ndat$old_label <- Ndat$tip_label
Ndat$tip_label <- Ndat$new_label
Ntree.root$tip.label <- Ndat$tip_label


colz = c("Sarbecovirus" = "darkorchid1", "Embecovirus"="darkgoldenrod1", "Gammacoronavirus" = "black", "Hibecovirus" = "royalblue", "Nobecovirus" = "tomato", "Merbecovirus" = "mediumseagreen")
shapez = c("bat host" =  24, "non-bat host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")



#first, get the legend.

pleg <- ggtree(Stree.root) %<+% Sdat + 
  geom_tippoint(aes(color=subgroup, shape=bat_host)) +
  geom_nodelab(size=1,nudge_x = -.08, nudge_y = .5) +
  scale_color_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlim(c(0,4))
pleg

legall <- cowplot::get_legend(pleg)


pS <- ggtree(Stree.root) %<+% Sdat + 
  geom_tippoint(aes(fill=subgroup, shape=bat_host), show.legend = F, size=3) +
  geom_nodelab(size=2.5,nudge_x = -.08, nudge_y = .5) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=3.4, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  geom_treescale(fontsize=4, x=.7,y=38, linesize = .5) + 
  #theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4.5))
pS


pE <- ggtree(Etree.root) %<+% Edat + 
  geom_tippoint(aes(fill=subgroup, shape=bat_host), show.legend = F, size=3) +
  geom_nodelab(size=2.5,nudge_x = -.08, nudge_y = .5) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=3.4, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  geom_treescale(fontsize=4, x=.9,y=38, linesize = .5) + 
  #theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4.5))
pE


pM <- ggtree(Mtree.root) %<+% Mdat + 
  geom_tippoint(aes(fill=subgroup, shape=bat_host), show.legend = F, size=3) +
  geom_nodelab(size=2.5,nudge_x = -.08, nudge_y = .5) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=3.4, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  geom_treescale(width=.2, fontsize=4, x=.5,y=38, linesize = .5) + 
  #theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4.5))
pM


pN <- ggtree(Ntree.root) %<+% Ndat + 
  geom_tippoint(aes(fill=subgroup, shape=bat_host), show.legend = F, size=3) +
  geom_nodelab(size=2.5,nudge_x = -.06, nudge_y = .7) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=3.4, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  geom_treescale(fontsize=4, x=.5,y=38, linesize = .5) + 
  #theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4.5))
pN


#and all together

pAminoAcid <- cowplot::plot_grid(pS, pE, pM, pN, nrow = 1, ncol = 4, labels = c("(A)", "(B)", "(C)", "(D)"), label_size = 20)

Fig4 <- cowplot::plot_grid(pAminoAcid, legall, nrow=2, ncol=1, rel_heights = c(1,.1))

ggsave(file = paste0(homewd, "/final-figures/Fig4.png"),
       plot = Fig4,
       units="mm",  
       width=180, 
       height=80, 
       #limitsize = F,
       scale=4)#, 
