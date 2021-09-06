rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd="/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig5/"))


#first, amino acid similarity


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

# 
# genome.df <- data.frame(position = c(1, 4316,
#                                      4316, 7002,
#                                      7002, 8324,
#                                      8324, 8580,
#                                      8580, 8659,
#                                      8659, 8881,
#                                      8881, 9378,
#                                      9378, 9820), 
#                         gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))


genome.df <- data.frame(position = c(1, 7002,
                                     7002, 8324,
                                     8324, 8580,
                                     8580, 8659,
                                     8659, 8881,
                                     8881, 9378,
                                     9378, 9820), 
                        gene = rep(c("ORF1ab", "S", "NS3", "E", "M", "N", "NS7"), each=2))


genome.df$gene <- factor(genome.df$gene, levels = unique(genome.df$gene))

#colz= c("HKU9"="firebrick3", "GCCDC1"="magenta", "GX2018.BtCoV92" = "purple", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")
colz= c("HKU9"="firebrick3", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")

p1leg <- ggplot(long.sim) + geom_line(aes(x=CenterPos, y=value, color=strain), size=1) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black") +
  facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz) + coord_cartesian(ylim=c(-.1,1)) +
  scale_x_continuous(breaks=c(0,10000/3.055,20000/3.055,30000/3.055), labels = c(0,10000, 20000,30000))

p1leg




leg1 <- cowplot::get_legend(p1leg)

p1 <- ggplot(long.sim) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.9) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black", show.legend = F) +
  facet_grid(~Query) + theme_bw() + xlab("genome position") + ylab("amino acid similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"),
        plot.margin = unit(c(.1, .1,0, .2), "cm"),
        legend.text = element_text(face="italic"),
        axis.text.y = element_text(size=12), axis.title.y = element_text(size=14),
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_color_manual(values=colz)
p1


simplot2 <- read.csv(file = "all_simplot.csv", header = T, stringsAsFactors = F)
head(simplot2)

#move to long
long.sim2 <- melt(simplot2, id.vars = c("CenterPos", "Alternate_ID", "Query"), measure.vars = c("HKU9", "Eidolon_helvum", "Alternate"))

head(long.sim2)
unique(long.sim2$variable)
long.sim2$variable <- as.character(long.sim2$variable)
long.sim2$variable[long.sim2$variable=="Alternate"] <- long.sim2$Alternate_ID[long.sim2$variable=="Alternate"] 
unique(long.sim2$variable)
names(long.sim2)[names(long.sim2)=="variable"] <- "strain"

long.sim2$strain[long.sim2$strain=="Eidolon_helvum"] <- "E. helvum bat coronavirus"
long.sim2$strain[long.sim2$strain=="Pteropus_rufus"] <- "P. rufus Nobecovirus"
long.sim2$strain[long.sim2$strain=="Rousettus_madagascariensis"] <- "R. madagascariensis Nobecovirus"

long.sim2$strain <- factor(long.sim2$strain, levels = c("HKU9", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))

#and plot

long.sim2$value[long.sim2$value<0] <- 0
long.sim2$Query[long.sim2$Query=="Pteropus_rufus"] <- "Pteropus rufus"
long.sim2$Query[long.sim2$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"

genome.df2 <- data.frame(position = c(284, 21400, 
                                     21311, 25447, 
                                     25447, 26189, 
                                     26189, 26422, 
                                     26422, 27114, 
                                     27186, 28738, 
                                     28828, 30263), 
                        gene = rep(c("Orf1ab", "S", "NS3", "E", "M", "N", "NS7"), each=2))


genome.df2$gene <- factor(genome.df2$gene, levels = unique(genome.df2$gene))



p2 <- ggplot(long.sim2) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.9) +
      geom_ribbon(data=genome.df2, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black", show.legend = F) +
      facet_grid(~Query) + theme_bw() + ylab("nucleotide similarity") +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text=element_blank(),
            plot.margin = unit(c(.1, .1,0, .2), "cm"),
            legend.text = element_text(face="italic"),
            axis.text.y = element_text(size=12), axis.title.y = element_text(size=14),
            axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_color_manual(values=colz)
p2


#add in the bootscan
bootplot <- read.csv(file = "all_bootscan.csv", header = T, stringsAsFactors = F)
head(bootplot)

#move to long
long.boot <- melt(bootplot, id.vars = c("CenterPos", "ID_Alternative", "Query"), measure.vars = c("HKU9", "Eidolon_helvum", "Alternative"))

head(long.boot)
unique(long.boot$variable)
long.boot$variable <- as.character(long.boot$variable)
long.boot$variable[long.boot$variable=="Alternative"] <- long.boot$ID_Alternative[long.boot$variable=="Alternative"] 
unique(long.boot$variable)
names(long.boot)[names(long.boot)=="variable"] <- "strain"

long.boot$strain_label <- NA

#ylab(bquote("r"^"*"~", virus growth")) 

long.boot$strain[long.boot$strain=="Eidolon_helvum"] <- "E. helvum bat coronavirus"
long.boot$strain[long.boot$strain=="Pteropus_rufus"] <- "P. rufus Nobecovirus"
long.boot$strain[long.boot$strain=="Rousettus_madagascariensis"] <-  "R. madagascariensis Nobecovirus"

long.boot$strain <- factor(long.boot$strain, levels = c("HKU9", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))

long.boot$value[long.boot$value<0] <- 0
long.boot$Query[long.boot$Query=="Pteropus_rufus"] <- "Pteropus rufus"
long.boot$Query[long.boot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"



p3 <- ggplot(long.boot) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.9) +
  geom_ribbon(data=genome.df2, aes(x=position, ymin=-10, ymax=-5, fill=gene), color="black", show.legend = F) +
  facet_grid(~Query) + theme_bw() + xlab("genome position") + ylab("% of Permuted Trees") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = unit(c(.1, .1,.1, .3), "cm"),
        legend.text = element_text(face="italic"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz)

p3







#and together
Fig5top <- cowplot::plot_grid(p1,p2,p3, nrow=3, ncol = 1, labels= c("(A)", "(B)", "(C)"), label_size = 16, label_x = -.01, rel_heights = c(1,.9,1.2))

#Fig5 <- cowplot::plot_grid(Fig5top, leg1, nrow = 1, ncol = 2, rel_widths = c(1,.2))


Fig5 <- cowplot::plot_grid(Fig5top, leg1, nrow = 2, ncol = 1, rel_heights = c(1,.1))

ggsave(file = paste0(homewd, "/final-figures/Fig5.png"),
       plot=Fig5,
       units="mm",  
       width=95, 
       height=70, 
       #limitsize = F,
       scale=4)#, 
