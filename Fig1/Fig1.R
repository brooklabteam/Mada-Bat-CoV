rm(list=ls())

#packages
library(sf)
library(mapplots)
library(scatterpie)
library(maptools)
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(ggrepel)

# To run this script, change the "mainwd" to wherever this folder
# ("Mada-Bat-CoV") is stored on your computer
# Also, make sure to download/clone the "Mada-GIS" folder to 
# your home computer. I recommend putting it in the same parent 
# directory as "Mada-Bat-CoV"

# For example, my two folders are stored at:

# "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/     ...AND
# "/Users/caraebrook/Documents/R/R_repositories/Mada-GIS/

# I keep all my github repos under "R_repositories"

#####################################################################
#####################################################################
# Set wd to data on this computer. Also ID homewd, assuming that 
# Mada-GIS is cloned to the same series of sub-folders
homewd = "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV/" 
#should be wherever "Mada-Bat-CoV" is stored on your home computer
basewd = paste(strsplit(homewd, "/")[[1]][1:6], collapse = "/")
mapwd = paste0(basewd, "/", "Mada-GIS")
setwd(paste0(homewd, "/", "Fig1/"))



#import madagascar shapfile
name<- paste0(mapwd, "/", "MDG-3/MDG_adm3.shp")
otl_file <- paste(name, sep="") 
orotl_shp <- st_read(otl_file)
#View(orotl_shp)  # Open attribute table
class(orotl_shp)

###import and configuration
# plot mada
# note that this may bog your computer down : I only 
# recommend printing it once to check. If too slow, you can always
# comment out the "print" line and save it temporarily as a pdf instead
# (save script is commented out below the plot)

p1<-ggplot() +  
  geom_sf(color = "lightgoldenrod1", fill = "lightgoldenrod1",data = orotl_shp)+
  coord_sf(xlim = c(42, 60), ylim = c(-26, -11.5), expand = FALSE)+
  theme_bw()+
  xlab("Longitude") + ylab("Latitude") 
#print(p1)
# 
#  ggsave(file = "tmp_map_1.pdf",
#         plot = p1,
#         units="mm",  
#         width=40, 
#         height=60, 
#         scale=3, 
#         dpi=300)
# 


#import CoV data
dat <- read.csv(file = paste0(homewd,"/metadata/all_NGS_8_3_2021_distribute.csv"), header = T, stringsAsFactors = F )
head(dat)
names(dat)

#add age class
#clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat$young_of_year)
dat$age_class <- dat$bat_age_class
dat$age_class[dat$age_class=="P" | dat$age_class=="L"] <- "A"
dat$age_class[dat$age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$age_class[dat$young_of_year=="yes"] <- "J"

# now subset the data to just include the columns of interest

dat <- dplyr::select(dat,roost_site,latitude_s, longitude_e,
                       collection_date, age_class, bat_sex,
                       species, sampleid, CoV)

head(dat)
unique(dat$roost_site)

#get sites
coordinate <- ddply(dat, .(roost_site), summarise, latitude_s=unique(latitude_s), longitude_e=unique(longitude_e))
coordinate <-subset(coordinate, roost_site=="Ambakoana" | roost_site=="AngavoKely" | roost_site=="Maromizaha")
coordinate$species <- c("Pteropus rufus", "Eidolon dupreanum", "Rousettus madagascariensis")
head(coordinate)

#plot sites on map
p2<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="#97B5CC",size=1,data=dat)+
  annotation_scale(location = "bl", width_hint = 0.05) +    # scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.02, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)
#print(p2)

#  ggsave(file = "tmp_map_2.pdf",
#         plot = p2,
#          units="mm",  
#          width=40, 
#          height=60, 
#          scale=3, 
#          dpi=300)
# # 

#load GPS point and label
p2b<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=species),
            fontface="italic",
            color = "#1B262C", size=3,
            nudge_x = c(-1,-2.5,4),
            nudge_y = c(1.5,.5,0),
            check_overlap = T)+
  annotation_scale(location = "bl", width_hint = 0.05) +    #scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.03, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)+
  geom_text_repel(segment.colour="black")+
  theme_bw() +theme(panel.grid = element_blank(), 
                    plot.title = element_text(color="black", size=12, face="bold"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.26,.90),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.background = element_rect(color="gray",size = .1),
                    legend.text = element_text(size = 9,face = "italic"))
#print(p2b)
# # # 
#     ggsave(file = "tmp_map_2b.pdf",
#            plot = p2b,
#            units="mm",  
#            width=40, 
#            height=60, 
#            scale=3, 
#            dpi=300)
# # # 


#plot one site per species
dat$roost_site[dat$species=="Pteropus rufus"] <- "Ambakoana"
dat$roost_site[dat$species=="Eidolon dupreanum"] <- "AngavoKely"
dat$roost_site[dat$species=="Rousettus madagascariensis"] <- "Maromizaha"

dat$longitude_e[dat$roost_site=="Ambakoana"] <- coordinate$longitude_e[coordinate$roost_site=="Ambakoana"]
dat$longitude_e[dat$roost_site=="AngavoKely"] <- coordinate$longitude_e[coordinate$roost_site=="AngavoKely"]
dat$longitude_e[dat$roost_site=="Maromizaha"] <- coordinate$longitude_e[coordinate$roost_site=="Maromizaha"]


dat$latitude_s[dat$roost_site=="Ambakoana"] <- coordinate$latitude_s[coordinate$roost_site=="Ambakoana"]
dat$latitude_s[dat$roost_site=="AngavoKely"] <- coordinate$latitude_s[coordinate$roost_site=="AngavoKely"]
dat$latitude_s[dat$roost_site=="Maromizaha"] <- coordinate$latitude_s[coordinate$roost_site=="Maromizaha"]




###Grouping data for scatterpie
dat$plot_class <- NA
dat$plot_class[dat$age_class=="J" & dat$CoV==1] <- "juvenile: CoV pos"
dat$plot_class[dat$age_class=="J" & dat$CoV==0] <- "juvenile: CoV neg"
dat$plot_class[dat$age_class=="A" & dat$CoV==1] <- "adult: CoV pos"
dat$plot_class[dat$age_class=="A" & dat$CoV==0] <- "adult: CoV neg"

pies <- ddply(dat, .(species, roost_site, latitude_s, longitude_e, plot_class), summarise, value=length(sampleid))

tot_sum = ddply(pies,.(species), summarise,N=sum(value))

pies <- merge(pies, tot_sum, by="species", all.x=T)

pies$plot_class <- factor(pies$plot_class, levels=c( "juvenile: CoV neg", "adult: CoV neg", "juvenile: CoV pos", "adult: CoV pos"))
###Get the pie data in the right format###
colz = c('adult: CoV neg' ="dodgerblue4", 'adult: CoV pos' ="firebrick4", 'juvenile: CoV neg' ="dodgerblue", 'juvenile: CoV pos' ="firebrick1")


p3<-ggplot() + 
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  scale_fill_manual(values=colz)



#print(p3)
# # 
#    ggsave(file = "tmp_map_3.pdf",
#           plot = p3,
#           units="mm",  
#           width=40, 
#           height=60, 
#           scale=3, 
#           dpi=300)
# #  
# 
# # copie of latitude (x.) and longitude (y.)
 pies$x2 <- pies$longitude_e
 pies$y2 <- pies$latitude_s
# 
# #manually move the pie chart in case there is an overlap (change x and y)
# 
 pies$x2[pies$species== "Pteropus rufus"] <- pies$longitude_e[pies$species== "Pteropus rufus"] + 1.5
 pies$y2[pies$species== "Pteropus rufus"] <- pies$latitude_s[pies$species== "Pteropus rufus"] + 1.5
 
 
 pies$x2[pies$species== "Eidolon dupreanum"] <- pies$longitude_e[pies$species== "Eidolon dupreanum"] - 1.5
 pies$y2[pies$species== "Eidolon dupreanum"] <- pies$latitude_s[pies$species== "Eidolon dupreanum"] - 1.5
 
 pies$x2[pies$species== "Rousettus madagascariensis"] <- pies$longitude_e[pies$species== "Rousettus madagascariensis"] + 1.5
 pies$y2[pies$species== "Rousettus madagascariensis"] <- pies$latitude_s[pies$species== "Rousettus madagascariensis"] - 1.5
 
 head(pies)

#plot pie chart 
#loko<-c("Rousettus madagascariensis"="#B200ED","Eidolon dupreanum"="#7FFF00","Pteropus rufus"="#0000FF")

p4 <- p2b+annotate("segment", x=pies$longitude_e, xend=pies$x2,y=pies$latitude_s,yend=pies$y2,size=.7,alpha=.5)+ # put the lines
  geom_scatterpie(aes(x=x2, y=y2, r=(N/100)), 
                          data = pies, cols="plot_class", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(0,.6,0,.3),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.77,.77),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 7.5)) +
  scale_fill_manual(values=colz) +
  geom_scatterpie_legend(c(50,100,150)/100,
                         x=54, y=-22.5, 
                         n=3,
                         labeller = function(x) paste(x*100,"indiv"))

#print(p4)

ggsave(file = paste0(homewd, "final-figures/Fig1.pdf"),
       plot=p4,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)

