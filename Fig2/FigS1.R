rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)

homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/"

#load txt file with coverage for each genome
Pruf <- read.delim(file = paste0(homewd, "Fig2/P_rufus.coverage"), header=F)
head(Pruf)
names(Pruf) <- c("genome", "position", "read_depth")

#and the others
R_mad240 <- read.delim(file = paste0(homewd, "Fig2/R_mada_MIZ240.coverage"), header=F)
head(R_mad240)
names(R_mad240) <- c("genome", "position", "read_depth")


R_mad178 <- read.delim(file = paste0(homewd, "Fig2/R_mada_MIZ178.coverage"), header=F)
head(R_mad178)
names(R_mad178) <- c("genome", "position", "read_depth")

#and separate out
R_mad178$sequence <- "Rousettus madagascariensis MIZ178"
R_mad240$sequence <- "Rousettus madagascariensis MIZ240"
Pruf$sequence <- "Pteropus rufus AMB130"

all.seq <- rbind(Pruf, R_mad178, R_mad240)
all.seq$sequence <- factor(all.seq$sequence, levels = c("Pteropus rufus AMB130", "Rousettus madagascariensis MIZ178", "Rousettus madagascariensis MIZ240"))
head(all.seq)

all.seq$label <- as.character(all.seq$sequence)
all.seq$label[all.seq$label=="Pteropus rufus AMB130"] <- "(A)"
all.seq$label[all.seq$label=="Rousettus madagascariensis MIZ178"] <- "(B)"
all.seq$label[all.seq$label=="Rousettus madagascariensis MIZ240"] <- "(C)"

p1 <- ggplot(data=all.seq) + 
      geom_line(aes(x=position, y=read_depth)) + 
      geom_label(aes(x=1,y=1, label=label)) +
      geom_ribbon(aes(x=position, ymin=0, ymax=read_depth), fill="gray") + 
      theme_bw() +
      facet_grid(sequence~.) + scale_y_log10() +
      ylab("read depth") +
      xlab("genome position") +
      theme(panel.grid = element_blank(),
            strip.text = element_text(face="italic", size=14),
            strip.background = element_rect(fill="white"),
            plot.margin = unit(c(.1, .1,.1, .1), "cm"),
            legend.text = element_text(face="italic"),
            axis.text = element_text(size=12), axis.title = element_text(size=14))
            

ggsave(file = paste0(homewd, "/final-figures/FigS1.png"),
       plot=p1,
       units="mm",  
       width=95, 
       height=70, 
       #limitsize = F,
       scale=4)#, 
