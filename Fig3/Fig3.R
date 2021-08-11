rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"

setwd(paste0(homewd, "/Fig3"))

#load the fig3a tree
treeA <-  read.tree(file = paste0(homewd, "Fig3/A-full-genome-phylogeny/Fig3A-raxml-output/Fig3A.raxml.supportTBE"))

#root it


rooted.tree.A <- root(treeA, which(treeA$tip.label == "NC_010800_1_|Turkey_coronavirus__complete_genome"))
#take a quick look in base R
plot(rooted.tree.A)

#load tree data prepared from elsewhere
dat <- read.csv(file = paste0(homewd, "Fig3/A-full-genome-phylogeny/fig3a_metadata_manual.csv"), header = T, stringsAsFactors = F)
head(dat)

#check subgroup names
unique(dat$sub_group)

colz = c("Alphacoronavirus" = "magenta", "Sarbecovirus" = "darkorchid1", "unclassified-Betacoronavirus"= "magenta", "Embecovirus"="darkgoldenrod1", "Gammacoronavirus" = "black", "Hibecovirus" = "royalblue", "Nobecovirus" = "tomato", "New-Madagascar-Nobecovirus" = "red", "Merbecovirus" = "mediumseagreen")

#pick order for the labels
dat$sub_group <- factor(dat$sub_group, levels = c("Alphacoronavirus", "Embecovirus", "Sarbecovirus", "unclassified-Betacoronavirus", "Hibecovirus", "Nobecovirus", "New-Madagascar-Nobecovirus", "Merbecovirus", "Gammacoronavirus"))   

#and add a "novel" category
dat$novel = 0
dat$novel[dat$sub_group=="New-Madagascar-Nobecovirus"] <- 1
dat$novel <- as.factor(dat$novel)

rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)

#take a glance
p <- ggtree(rooted.tree.A) %<+% dat + geom_tippoint(aes(fill=sub_group)) +
      geom_tiplab(size=1) + geom_nodelab(size=1) +
    scale_fill_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$tip_label <- paste(dat$accession_num, "|", dat$host, "|", dat$country, "|", dat$collection_year)
#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.dat <- data.frame(old_tip_label=rooted.tree.A$tip.label, num =1:length(rooted.tree.A$tip.label))
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)
rooted.tree.A$tip.label <- tree.dat$tip_label

dat$bat_host[dat$bat_host==0] <- "non-bat host"
dat$bat_host[dat$bat_host==1] <- "bat host"
dat$bat_host <- as.factor(dat$bat_host)
shapez = c("bat host" =  24, "non-bat host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")


p1 <- ggtree(rooted.tree.A) %<+% dat + geom_tippoint(aes(fill=sub_group, shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4))
p1

# # #save temporary draft, now make figure 3B
#  ggsave(file = paste0(homewd, "/Fig3/A-full-genome-phylogeny/Fig3A-draft.png"),
#         units="mm",  
#         width=100, 
#         height=100, 
#         #limitsize = F,
#         scale=4)#, 


#now add in details for the second panel

#### And Fig 3 B
#load the fig3b tree
treeB <-  read.tree(file = paste0(homewd, "Fig3/B-RdRP-phylogeny/Fig3B-raxml-output/Fig3B.raxml.supportTBE"))


#take a quick look in base R
plot(treeB)


#root the tree in the outgroup
rooted.tree.B <- root(treeB, which(treeB$tip.label == "GammaCoV_NC_010800_1_Turkey"))
p3Broot <- ggtree(rooted.tree.B) + geom_tiplab() # rooted by the outgroup

p3Broot


#load tree data prepared from elsewhere
datB <- read.csv(file = paste0(homewd, "Fig3/B-RdRP-phylogeny/fig3b_metadata_manual.csv"), header = T, stringsAsFactors = F)
head(datB)

#check subgroup names
unique(datB$sub_group)

colz = c("Alphacoronavirus" = "magenta", "Sarbecovirus" = "darkorchid1", "unclassified-Betacoronavirus"= "magenta", "Embecovirus"="darkgoldenrod1", "Gammacoronavirus" = "black", "Hibecovirus" = "royalblue", "Nobecovirus" = "tomato", "New-Madagascar-Nobecovirus" = "red", "Merbecovirus" = "mediumseagreen")

#pick order for the labels
datB$sub_group <- factor(datB$sub_group, levels = c("Alphacoronavirus", "Embecovirus", "Sarbecovirus", "unclassified-Betacoronavirus", "Hibecovirus", "Nobecovirus", "New-Madagascar-Nobecovirus", "Merbecovirus", "Gammacoronavirus"))   


rooted.tree.B$node.label <- round(as.numeric(rooted.tree.B$node.label)*100, 0)

#take a glance in the same manner as in A
pB <- ggtree(rooted.tree.B) %<+% datB + geom_tippoint(aes(fill=sub_group)) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz) + theme(legend.position = c(.2,.85), legend.title = element_blank())
pB #looks great

#now get new tip labels
datB$old_tip_label <- datB$tip_label
datB$tip_label <- paste(datB$accession_num, "|", datB$host, "|", datB$country, "|", datB$collection_year)
#after Gwen checks these, can later manually edit any that are "NA"

#make sure to sort in order
tree.datB <- data.frame(old_tip_label=rooted.tree.B$tip.label, num =1:length(rooted.tree.B$tip.label))
tree.datB <- merge(tree.datB, datB, by = "old_tip_label", all.x = T, sort = F)
rooted.tree.B$tip.label <- tree.datB$tip_label

datB$bat_host[datB$bat_host==0] <- "non-bat host"
datB$bat_host[datB$bat_host==1] <- "bat host"
datB$bat_host <- as.factor(datB$bat_host)

datB$novel = 0
datB$novel[datB$sub_group=="New-Madagascar-Nobecovirus"] <- 1
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

#save draft
# 
# ggsave(file = paste0(homewd, "/Fig3/B-RdRP-phylogeny/Fig3B-draft.png"),
#        units="mm",  
#        width=50, 
#        height=50, 
#        #limitsize = F,
#        scale=3)#, 
# 


#####now combine the two together somehow


###wrking on p1

p1 <- ggtree(rooted.tree.A) %<+% dat + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), size=3) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=3, hjust = -.15) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(c(0,4))
p1

#collapse the alpha clade (all bat CoVs)
alpha_node = MRCA(rooted.tree.A, which(rooted.tree.A$tip.label == "NC_048211 | Suncus_murinus | China | 2015" ),which(rooted.tree.A$tip.label == "NC_018871 | bat | China | 2021" ))
p1.2.leg <- p1 %>% ggtree::collapse(node = alpha_node) +
  geom_tippoint(aes(color=sub_group, shape=bat_host), size=2) +
  scale_color_manual(values=colz) + 
  geom_nodelab(size=1,nudge_x = -.04, nudge_y = .7) +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.title = element_blank(), legend.text = element_text(size=12)) +
  geom_point2(aes(subset=(node==alpha_node)), shape=24, size=5, fill="magenta") 

p1.2.leg
#separate legend
leg.all <- cowplot::get_legend(p1.2.leg)

p1.3 <- p1 %>% ggtree::collapse(node = alpha_node) +
  geom_nodelab(size=2.5,nudge_x = -.04, nudge_y = .7) +
  theme(legend.position = "none", legend.title = element_blank()) +
  geom_point2(aes(subset=(node==alpha_node)), shape=24, size=7, fill="magenta") +
  xlim(c(0,3))

p1.3

#new p2
p2.1 <- ggtree(rooted.tree.B) %<+% datB + 
  geom_tippoint(aes(fill=sub_group, shape=bat_host), show.legend = F, size=5) +
  geom_nodelab(size=3,nudge_x = -.03, nudge_y = .4) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill()+
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=4, hjust = -.13) + 
  scale_fill_manual(values=colz2) + 
  theme(legend.position = "right", legend.title = element_blank()) +
  xlim(c(0,1.5))
p2.1 

#great, now need to flip some of the clases to match plot on the left

node_flip_embeco_Nobeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_048217 | Mus_musculus | NA | NA" ),which(rooted.tree.B$tip.label == "KP696742 | Pteropus_rufus | Madagascar | 2010" ))
node_flip_Merbeco_Nobeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_039207 | Erinaceus_europaeus | Germany | 2012" ),which(rooted.tree.B$tip.label == "KP696742 | Pteropus_rufus | Madagascar | 2010" ))
node_flip_Sarbeco_Nobeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_004718 | Homo_sapiens | Canada | NA" ),which(rooted.tree.B$tip.label == "KP696742 | Pteropus_rufus | Madagascar | 2010" ))
node_flip_Sarbeco_Hibeco = MRCA(rooted.tree.B, which(rooted.tree.B$tip.label == "NC_004718 | Homo_sapiens | Canada | NA" ),which(rooted.tree.B$tip.label == "NC_025217 | Hipposideros_pratti | China | 2013" ))

p2.2 <- p2.1 %>% ggtree::rotate(node = node_flip_embeco_Nobeco)
p2.3 <- p2.2 %>% ggtree::rotate(node = node_flip_Merbeco_Nobeco)
p2.4 <- p2.3 %>% ggtree::rotate(node = node_flip_Sarbeco_Nobeco)
p2.5 <- p2.4 %>% ggtree::rotate(node = node_flip_Sarbeco_Hibeco)

Fig3 <- cowplot::plot_grid(p1.3,p2.5, ncol=2, nrow=1, labels = c("A.", "B."), label_size = 22, label_x = .03, label_y = .98)

Fig3all <- cowplot::plot_grid(Fig3,leg.all, ncol=1, nrow=2, rel_heights = c(1,.1))



 ggsave(file = paste0(homewd, "/Fig3/Fig3-draft.png"),
        units="mm",  
        width=150, 
        height=100, 
        #limitsize = F,
        scale=4)#, 
