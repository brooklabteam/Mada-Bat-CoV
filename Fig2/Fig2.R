rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig2/"))

simplot <- read.csv(file = paste0(homewd, "/Fig2/all_simplotAA_noInsertion.csv"), header = T, stringsAsFactors = F)
head(simplot)



#move to long
#long.sim <- melt(simplot, id.vars = c("CenterPos", "Alternate_ID", "Query"), measure.vars = c("HKU9", "Eidolon_helvum", "GX2018.BtCoV92", "GCCDC1", "Alternate"))
long.sim <- melt(simplot, id.vars = c("CenterPos", "Alternate_ID", "Query"), measure.vars = c("HKU9", "Eidolon_helvum", "Alternate"))


head(long.sim)
unique(long.sim$variable)
long.sim$variable <- as.character(long.sim$variable)
long.sim$variable[long.sim$variable=="Alternate"] <- long.sim$Alternate_ID[long.sim$variable=="Alternate"] 
unique(long.sim$variable)
names(long.sim)[names(long.sim)=="variable"] <- "strain"

long.sim$strain[long.sim$strain=="Eidolon_helvum"] <- "E. helvum bat coronavirus"
long.sim$strain[long.sim$strain=="Pteropus_rufus"] <- "P. rufus Nobecovirus"
long.sim$strain[long.sim$strain=="R_madagascariensis"] <- "R. madagascariensis Nobecovirus"

#long.sim$strain <- factor(long.sim$strain, levels = c("HKU9", "GCCDC1", "GX2018.BtCoV92", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))
long.sim$strain <- factor(long.sim$strain, levels = c("HKU9", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))


#and plot

long.sim$value[long.sim$value<0] <- 0
long.sim$Query[long.sim$Query=="Pteropus_rufus"] <- "Pteropus rufus"
long.sim$Query[long.sim$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"
long.sim$Query <- factor(long.sim$Query, labels = c("Pteropus rufus", "Rousettus madagascariensis"))

# genome.df <- data.frame(position = c(1, 4390,
#                                      4390, 7064,
#                                      7064, 8417,
#                                      8417, 8654,
#                                      8654, 8740,
#                                      8740, 8959,
#                                      8959, 9452,
#                                      9453, 10222), 
#                         gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))

genome.df <- data.frame(position = c(1, 4316,
                                     4316, 7002,
                                     7002, 8324,
                                     8324, 8580,
                                     8580, 8659,
                                     8659, 8881,
                                     8881, 9378,
                                     9378, 9820), 
                        gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))



genome.df$gene <- factor(genome.df$gene, levels = unique(genome.df$gene))

#colz= c("HKU9"="firebrick3", "GCCDC1"="magenta", "GX2018.BtCoV92" = "purple", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")
colz= c("HKU9"="firebrick3", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")

p1leg <- ggplot(long.sim) + geom_line(aes(x=CenterPos, y=value, color=strain), size=1) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black") +
  facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"), 
        legend.position = "right", #legend.direction = "horizontal",
        legend.text = element_text(face="italic"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz) + coord_cartesian(ylim=c(-.1,1)) +
  scale_x_continuous(breaks=c(0,10000/3.055,20000/3.055,30000/3.055), labels = c(0,10000, 20000,30000))

p1leg




leg1 <- cowplot::get_legend(p1leg)

p1 <- ggplot(long.sim) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.7) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black", show.legend = F) +
  facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(face="italic"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz)
p1




p2 <- ggplot(long.sim) + geom_bar(aes(x=CenterPos, y=value, color=strain), stat = "identity", position = "nudge", size=.7) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black", show.legend = F) +
  facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"),
        legend.text = element_text(face="italic"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz)
p2
