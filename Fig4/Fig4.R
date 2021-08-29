rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd="/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig4/"))

simplot <- read.csv(file = "all_simplot.csv", header = T, stringsAsFactors = F)
head(simplot)

#move to long
long.sim <- melt(simplot, id.vars = c("CenterPos", "Alternate_ID", "Query"), measure.vars = c("HKU9", "Eidolon_helvum", "Alternate"))

head(long.sim)
unique(long.sim$variable)
long.sim$variable <- as.character(long.sim$variable)
long.sim$variable[long.sim$variable=="Alternate"] <- long.sim$Alternate_ID[long.sim$variable=="Alternate"] 
unique(long.sim$variable)
names(long.sim)[names(long.sim)=="variable"] <- "strain"

long.sim$strain[long.sim$strain=="Eidolon_helvum"] <- "Eidolon helvum"
long.sim$strain[long.sim$strain=="Pteropus_rufus"] <- "Pteropus rufus"
long.sim$strain[long.sim$strain=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"

long.sim$strain <- factor(long.sim$strain, levels = c("HKU9", "Eidolon helvum", "Pteropus rufus", "Rousettus madagascariensis"))

#and plot

long.sim$value[long.sim$value<0] <- 0
long.sim$Query[long.sim$Query=="Pteropus_rufus"] <- "Pteropus rufus"
long.sim$Query[long.sim$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"

genome.df <- data.frame(position = c(284, 21400, 
                                     21311, 25447, 
                                     25447, 26189, 
                                     26189, 26422, 
                                     26422, 27114, 
                                     27186, 28738, 
                                     28828, 30263), 
                        gene = rep(c("Orf1ab", "S", "NS3", "E", "M", "N", "NS7"), each=2))

genome.df$gene <- factor(genome.df$gene, levels = unique(genome.df$gene))
colz= c("HKU9"="firebrick3", "Eidolon helvum" = "royalblue", "Pteropus rufus" = "purple3", "Rousettus madagascariensis" = "forestgreen")
p1leg <- ggplot(long.sim) + geom_line(aes(x=CenterPos, y=value, color=strain), size=1) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black") +
  facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("similarity") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"), 
        legend.position = "right", #legend.direction = "horizontal",
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz)
p1leg

leg1 <- cowplot::get_legend(p1leg)

p1 <- ggplot(long.sim) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.7) +
      geom_ribbon(data=genome.df, aes(x=position, ymin=-.1, ymax=-.05, fill=gene), color="black", show.legend = F) +
      facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("similarity") +
      theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
            strip.background = element_rect(fill="white"),
            axis.text = element_text(size=12), axis.title = element_text(size=14)) +
      scale_color_manual(values=colz)
p1


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

long.boot$strain[long.boot$strain=="Eidolon_helvum"] <- "Eidolon helvum"
long.boot$strain[long.boot$strain=="Pteropus_rufus"] <- "Pteropus rufus"
long.boot$strain[long.boot$strain=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"

long.boot$strain <- factor(long.boot$strain, levels = c("HKU9", "Eidolon helvum", "Pteropus rufus", "Rousettus madagascariensis"))

long.boot$value[long.boot$value<0] <- 0
long.boot$Query[long.boot$Query=="Pteropus_rufus"] <- "Pteropus rufus"
long.boot$Query[long.boot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"



p2 <- ggplot(long.boot) + geom_line(aes(x=CenterPos, y=value, color=strain), show.legend = F, size=.7) +
  geom_ribbon(data=genome.df, aes(x=position, ymin=-10, ymax=-5, fill=gene), color="black", show.legend = F) +
  facet_grid(Query~.) + theme_bw() + xlab("genome position") + ylab("% of Permuted Trees") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14),
        strip.background = element_rect(fill="white"), 
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  scale_color_manual(values=colz)
p2


#and together
Fig4left <- cowplot::plot_grid(p1,p2, nrow=1, ncol = 2, labels= c("A.", "B."), label_size = 18, label_x = -.02)

Fig4 <- cowplot::plot_grid(Fig4left, leg1, nrow = 1, ncol = 2, rel_widths = c(1,.2))


ggsave(file = paste0(homewd, "/final-figures/Fig4.png"),
       units="mm",  
       width=90, 
       height=40, 
       #limitsize = F,
       scale=4)#, 
