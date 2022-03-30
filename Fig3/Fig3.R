rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

homewd= "/Users/carabrook/Developer/Mada-Bat-CoV/"

setwd(paste0(homewd, "/Fig3"))

#load the fig3a tree
treeA <-  read.tree(file = paste0(homewd, "Fig3/A-full-genome-phylogeny/3-Fig3A-raxml-output/Fig3A.raxml.supportFBP"))

#root it


rooted.tree.A <- root(treeA, which(treeA$tip.label == "NC_010800_1_Turkey_coronavirus"))
#take a quick look in base R
plot(rooted.tree.A)

#load tree data prepared from elsewhere
dat <- read.csv(file=paste0(homewd,"Fig3/A-full-genome-phylogeny/fig3a_allbetacov_metadata_manual.csv"), header = T, stringsAsFactors = F)
head(dat)
#check subgroup names
unique(dat$sub_group)

colz = c("Sarbecovirus" = "darkorchid1", "unclassified-Betacoronavirus"= "magenta", "Embecovirus"="darkgoldenrod1", "Gammacoronavirus" = "black", "Hibecovirus" = "royalblue", "Nobecovirus" = "tomato", "Merbecovirus" = "mediumseagreen")

#pick order for the labels
dat$sub_group <- factor(dat$sub_group, levels = c("Embecovirus", "Sarbecovirus", "unclassified-Betacoronavirus", "Hibecovirus", "Nobecovirus", "Merbecovirus", "Gammacoronavirus"))   

#and add a "novel" category
dat$novel = 0
dat$novel[dat$country=="Madagascar"] <- 1
dat$novel <- as.factor(dat$novel)

#rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p <- ggtree(rooted.tree.A) %<+% dat + geom_tippoint(aes(color=sub_group)) +
      geom_tiplab(size=1) + geom_nodelab(size=1) +
    scale_color_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA
dat$new_label[!is.na(dat$strain)] <- paste(dat$accession_num[!is.na(dat$strain)], " | ", 
                                           dat$strain[!is.na(dat$strain)], " | ", 
                                           dat$host[!is.na(dat$strain)], " | ",
                                           dat$country[!is.na(dat$strain)], " | ",
                                           dat$collection_year[!is.na(dat$strain)])

dat$new_label[is.na(dat$strain)] <- paste(dat$accession_num[is.na(dat$strain)], " | ", 
                                           dat$host[is.na(dat$strain)], " | ",
                                           dat$country[is.na(dat$strain)], " | ",
                                           dat$collection_year[is.na(dat$strain)])

#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.dat <- data.frame(old_tip_label=rooted.tree.A$tip.label, num =1:length(rooted.tree.A$tip.label))
head(tree.dat)
head(dat)
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, accession_num, strain, host, bat_host, country, collection_date, collection_year, sub_group, novel, old_tip_label)

rooted.tree.A$tip.label <- tree.dat$tip_label

tree.dat$bat_host[tree.dat$bat_host==0] <- "non-bat host"
tree.dat$bat_host[tree.dat$bat_host==1] <- "bat host"
tree.dat$bat_host <- as.factor(tree.dat$bat_host)
shapez = c("bat host" =  24, "non-bat host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")


p1 <- ggtree(rooted.tree.A) %<+% tree.dat + geom_tippoint(aes(fill=sub_group, shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4))
p1




#now add in details for the second panel

#### And Fig 3 B
#load the fig3b tree
treeB <-  read.tree(file = paste0(homewd, "Fig3/B-RdRP-phylogeny/3-Fig3B-raxml-output/Fig3B.raxml.supportFBP"))


#take a quick look in base R
plot(treeB)


#root the tree in the outgroup
rooted.tree.B <- root(treeB, which(treeB$tip.label == "GammaCoV_NC_010800_1_Turkey"))
p3Broot <- ggtree(rooted.tree.B) + geom_tiplab() # rooted by the outgroup


#load tree data prepared from elsewhere
datB <- read.csv(file = paste0(homewd, "Fig3/B-RdRP-phylogeny/fig3b_metadata_manual.csv"), header = T, stringsAsFactors = F)
head(datB)


#check subgroup names
unique(datB$sub_group)


colz = c("Sarbecovirus" = "darkorchid1", "unclassified-Betacoronavirus"= "magenta", "Embecovirus"="darkgoldenrod1", "Gammacoronavirus" = "black", "Hibecovirus" = "royalblue", "Nobecovirus" = "tomato", "Merbecovirus" = "mediumseagreen")


#pick order for the labels
datB$sub_group <- factor(datB$sub_group, levels = c("Embecovirus", "Sarbecovirus", "Hibecovirus", "Nobecovirus", "Merbecovirus", "Gammacoronavirus"))   


#rooted.tree.B$node.label <- round(as.numeric(rooted.tree.B$node.label)*100, 0)

#take a glance in the same manner as in A
pB <- ggtree(rooted.tree.B) %<+% datB + geom_tippoint(aes(fill=sub_group)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
pB #looks great

#now get new tip labels
datB$old_tip_label <- datB$tip_label
datB$tip_label <- NA
datB$tip_label[!is.na(datB$strain)] <-paste(datB$accession_num[!is.na(datB$strain)], " | ",
                                            datB$strain[!is.na(datB$strain)], " | ",
                                            datB$host[!is.na(datB$strain)], " | ", 
                                            datB$country[!is.na(datB$strain)], " | ",
                                            datB$collection_year[!is.na(datB$strain)])

datB$tip_label[is.na(datB$strain)] <-paste(datB$accession_num[is.na(datB$strain)], " | ",
                                            datB$host[is.na(datB$strain)], " | ", 
                                            datB$country[is.na(datB$strain)], " | ",
                                            datB$collection_year[is.na(datB$strain)])

#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.datB <- data.frame(old_tip_label=rooted.tree.B$tip.label, num =1:length(rooted.tree.B$tip.label))
tree.datB <- merge(tree.datB, datB, by = "old_tip_label", all.x = T, sort = F)
rooted.tree.B$tip.label <- tree.datB$tip_label

datB$bat_host[datB$bat_host==0] <- "non-bat host"
datB$bat_host[datB$bat_host==1] <- "bat host"
datB$bat_host <- as.factor(datB$bat_host)

datB$novel = 0
datB$novel[datB$accession_num=="OK020086" |
           datB$accession_num=="OK067321" |
           datB$accession_num=="OK067320" | 
           datB$accession_num=="OK020089" |
           datB$accession_num=="OK067319" |
           datB$accession_num=="OK020087" |
           datB$accession_num=="OK020088" ] <- 1

datB$novel <- as.factor(datB$novel)
colz2 = c('1' =  "yellow", '0' = "white")

shapez = c("bat host" =  24, "non-bat host" = 21)

p2 <- ggtree(rooted.tree.B) %<+% datB + geom_tippoint(aes(fill=sub_group, shape=bat_host)) +
  geom_nodelab(size=1.5,nudge_x = -.02, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) +
  new_scale_fill()+
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  xlim(c(0,2))
p2



#####now combine the two together somehow


###wrking on p1

p1 <- ggtree(rooted.tree.A) %<+% tree.dat + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), size=3, show.legend = F) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = .9) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.08) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=1,y=124, linesize = .5) + 
  xlim(c(0,4))
p1

#and flip some clades
node_flip_Embeco_Merbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.A$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))
node_flip_Merbeco_Sarbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.A$tip.label == "MZ081380  |  SARSr_CoV  |  Rhinolophus_stheno  |  China  |  2020"))
node_flip_Embeco_Sarbeco1 = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "MZ081380  |  SARSr_CoV  |  Rhinolophus_stheno  |  China  |  2020" ),which(rooted.tree.A$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))

p1.2 <- p1 %>% ggtree::rotate(node = node_flip_Embeco_Sarbeco1 )


#p1.2 <- p1 %>% ggtree::rotate(node = node_flip_Merbeco_Sarbeco1)
#p1.3 <- p1.2 %>% ggtree::rotate(node = node_flip_Embeco_Sarbeco1)


#collapse the alpha clade (all bat CoVs)
#alpha_node = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_048211 | Suncus_murinus | China | 2015" ),which(rooted.tree.A$tip.label == "NC_018871 | bat | China | 2021" ))

p1.2.leg <- ggtree(rooted.tree.A) %<+% tree.dat + 
  geom_tippoint(aes(color=sub_group, shape=bat_host), size=3) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = .9) +
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.08) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size=12)) +
  geom_treescale(fontsize=4, x=1,y=124, linesize = .5) + 
  xlim(c(0,4))
p1.2.leg

#separate legend
leg.all <- cowplot::get_legend(p1.2.leg)


#new p2
p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  xlim(c(0,1.5))
p2.1 

#add lineage clade labels bars

#nodebase
clade.a <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "EF065514  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005" ),which(rooted.tree.B$tip.label == "HM211098  |  HKU9  |  Rhinolophus_sinicus  |  China  |  2005" ))
clade.b <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "OK020089  |  Rousettus_madagascariensis  |  Madagascar  |  2018" ),which(rooted.tree.B$tip.label == "MG693172  |  Eidolon_helvum  |  Cameroon  |  2013" ))
clade.c <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_030886  |  GCCDC1  |  Rousettus_leschenaulti  |  China  |  2014" ),which(rooted.tree.B$tip.label == "MT350598  |  GCCDC1  |  Eonycteris_spelaea  |  Singapore  |  2016" ))
clade.d <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "MK211379  |  GX2018  |  Rhinolophus_affinis  |  China  |  2016" ),which(rooted.tree.B$tip.label == "MK492263  |  BatCoV92  |  Cynopteris_brachyotis  |  Singapore  |  2015" ))
clade.e <- MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "OK020087  |  Pteropus_rufus  |  Madagascar  |  2018" ),which(rooted.tree.B$tip.label == "OK067319  |  Pteropus_rufus  |  Madagascar  |  2018" ))




p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=4, hjust = -.1) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  geom_treescale(fontsize=4, x=.3,y=50, linesize = .5) + 
  geom_cladelabel(node = clade.a, label = "HKU9", offset = 1, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.b, label = "atop(African,italic(Eidolon))", offset = .8, fontsize = 6.5, color="tomato", parse=T) +
  geom_cladelabel(node = clade.c, label = "GCCDC1", offset = 1.05, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.d, label = "BtCoV92 /\nGX2018", offset = 1.03, fontsize = 6.5, color="tomato", offset.text = .01) +
  geom_cladelabel(node = clade.e, label = "atop(Madagascar,italic(Pteropus))", offset = .7, fontsize = 6.5, color="tomato" , parse = T) +
  xlim(c(0,1.8))
p2.1 


#great, now need to flip some of the clases to match plot on the left

node_flip_Embeco_Merbeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.B$tip.label == "NC_006213  |  HCoV_OC43  |  Homo_sapiens  |  USA  |  1960"  ))
node_flip_Sarbeco_Hibeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_025217  |  Hipposideros_pratti  |  China  |  2013" ),which(rooted.tree.B$tip.label == "NC_004718  |  SARS_CoV  |  Homo_sapiens  |  Canada  |  2003" ))
node_flip_Embeco_Nobeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_019843  |  MERS  |  Homo_sapiens  |  Saudi_Arabia  |  2012" ),which(rooted.tree.B$tip.label == "EF065516  |  HKU9  |  Rousettus_leschenaulti  |  China  |  2005"   ))


p2.2 <- p2.1 %>% ggtree::rotate(node = node_flip_Embeco_Merbeco)
p2.3 <- p2.2 %>% ggtree::rotate(node = node_flip_Sarbeco_Hibeco)
p2.4 <- p2.3 %>% ggtree::rotate(node = node_flip_Embeco_Nobeco)


Fig3 <- cowplot::plot_grid(p1.2,p2.4, ncol=2, nrow=1, labels = c("(A)", "(B)"), label_size = 22, label_x = .03, label_y = .98)

Fig3all <- cowplot::plot_grid(Fig3,leg.all, ncol=1, nrow=2, rel_heights = c(1,.1))


#and save to the final figures

 ggsave(file = paste0(homewd, "/final-figures/Fig3.png"),
        units="mm",  
        width=150, 
        height=100, 
        #limitsize = F,
        scale=4)#, 

 