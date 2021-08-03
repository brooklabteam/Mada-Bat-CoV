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
  geom_sf(color = "#EAF4F8", fill = "#EAF4F8",data = orotl_shp)+
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
dat <- read.csv(file = paste0(homewd,"all_NGS_8_3_2021_distribute.csv"), header = T, stringsAsFactors = F )
head(dat)
names(dat)


# now subset the data to just include the columns of interest

data1 <- dplyr::select(dat,roost_site,latitude_s, longitude_e, 
                       collection_date,
                       bat_species, sampleid)

#and just the three fruit bat species
unique(dat$bat_species)
dat <- subset(dat, bat_species!="Hipposideros commersoni" & bat_species!="Mormopterus jugularis" & bat_species!="Asio madagascariensis")
head(dat)

#group all the ankarana sites together (and Moramanga)
dat$roost_site[dat$roost_site=="Ankarana_Canyon" |dat$roost_site=="Ankarana_Cathedral" | dat$roost_site=="Ankarana_Chauves_Souris"] <- "Ankarana"
dat$roost_site[dat$roost_site=="AngavoBe"  |dat$roost_site=="AngavoKely" |dat$roost_site=="Angavokely"] <- "Moramanga"
dat$roost_site[dat$roost_site=="Mangarivotra"| dat$roost_site=="Marovitsika"| dat$roost_site=="Maromizaha"] <- "Moramanga"
dat$roost_site[dat$roost_site=="Marotsipohy" |  dat$roost_site=="Ambakoana"| dat$roost_site=="Mahialambo"| dat$roost_site=="Lakato"] <- "Moramanga"
head(dat)

unique(dat$roost_site)

###import GPS (latitude and longitude for Ankarana, Makira, Moramanga, Mahabo)
coordinate<-read.csv("siteGPS.csv")
head(coordinate)

#group Ankarana GPS together (and Moramanga)
coordinate$roost_site[coordinate$roost_site=="Ankarana_Cathedral"] <- "Ankarana"
coordinate = subset(coordinate, roost_site !="Ankarana_Canyon" & roost_site != "Ankarana_Chauves_Souris" )

coordinate$roost_site[coordinate$roost_site=="Lakato"] <- "Moramanga"
coordinate = subset(coordinate, roost_site !="AngavoKely" & roost_site != "Angavokely" & roost_site != "AngavoBe" & roost_site != "Marovitsika" & roost_site != "Maromizaha" & roost_site != "Ambakoana" & roost_site != "Mahialambo")

head(coordinate)

#load GPS point
p2<-p1+geom_point(aes(x=x,y=y),color="#97B5CC",size=1,data=coordinate)+
  annotation_scale(location = "bl", width_hint = 0.05) +    # scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.02, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)
#print(p2)

 # ggsave(file = "tmp_map_2.pdf",
 #        plot = p2,
 #        units="mm",  
 #        width=40, 
 #        height=60, 
 #        scale=3, 
 #        dpi=300)
# 

#load GPS point and label
p2b<-p1+geom_point(aes(x=x,y=y),color="#651441",size=1,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=x, y=y, label=roost_site), 
            color = "#1B262C", size=3.5,
            nudge_x = c(-1.5,-2,0,-1.4,-2,-1.5),
            nudge_y = c(0,0,.5,0,.3,0),
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
# # 
#   ggsave(file = "tmp_map_2b.pdf",
#          plot = p2b,
#          units="mm",  
#          width=40, 
#          height=60, 
#          scale=3, 
#          dpi=300)
# # 


###Grouping data###
binned_sites <- dat %>%
  dplyr::select(roost_site, bat_species) %>%
  group_by(roost_site, bat_species) %>%
  summarize(n = n())


###Adding identifiers for merge###
dat <- as.data.frame(binned_sites)
head(dat)

dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)],as.factor)
colnames(dat)[[1]] <- "roost_site"
colnames(dat)[[2]] <- "bat_species"
colnames(dat)[[3]] <- "n_indiv"
head(dat)
dat$nombre<-dat$n_indiv # DUplication 


############################################################################################
#sums of individuals per site
############################################################################################

#calcul of radius 

dat1<-dat%>% 
  group_by(roost_site) %>% 
  mutate(somme = sum(nombre))

#View(dat1)

############################################################################################
##Sum of bats captured for all sites
###########################################################################################
#Proportion 
rehetra<-c()

for (i in length(dat1$roost_site)) { 
  rehetra<-append(rehetra, sum(dat1$n_indiv))
  print(rehetra)
}

dat1$rehetra<-rehetra
head(dat1)


###Get the pie data in the right format###
pies <- dat %>% 
  group_by(roost_site) %>% 
  mutate(total = n(),
         radius  = length(unique(bat_species)))%>%
  group_by(roost_site, bat_species) %>% 
  summarise(prop = n_indiv/total,
            radius = mean(radius)) %>% 
  distinct()

#View(pies)


#pie calcul
pies$prop <- pies$prop
pies$value <- dat1$n_indiv
pies$nombre<-dat1$nombre
pies$nombre<-dat1$somme
pies$rehetra<-dat1$rehetra

#now give the pies xy coordinates by mergeing with the "coordinate" dataset
pies <- merge(pies, coordinate, by = "roost_site", all.x=T, sort =F)
names(pies)
#View(pies)

head(pies)


# Calculation of proportion and new radius
pies$proportion<-pies$value/pies$rehetra*100 #Proportion
#View(pies)
pies$propt<-pies$value/pies$nombre

pies$rayon<-pies$nombre/pies$rehetra

###Get the pie data in the right format###
p3<-ggplot(data=pies) + 
  geom_scatterpie(aes(x=x, y=y, group = roost_site, r = log10(rayon)), 
                  data = pies, cols = "bat_species", long_format=TRUE, color = NA)

#print(p3)
# # 
#   ggsave(file = "tmp_map_3.pdf",
#          plot = p3,
#          units="mm",  
#          width=40, 
#          height=60, 
#          scale=3, 
#          dpi=300)
#  

# copie of latitude (x.) and longitude (y.)
pies$x. <- pies$x
pies$y. <- pies$y

#manually move the pie chart in case there is an overlap (change x and y)

pies$x. <- ifelse(pies$roost_site == "Moramanga", pies$x. + 3.2, pies$x.)
pies$y. <- ifelse(pies$roost_site == "Moramanga", pies$y. + 1, pies$y.)

pies$x. <- ifelse(pies$roost_site == "Mahabo", pies$x. -0.5, pies$x.)
pies$y. <- ifelse(pies$roost_site == "Mahabo", pies$y. -1, pies$y.)

pies$x. <- ifelse(pies$roost_site == "Makira", pies$x. - 0.7, pies$x.)
pies$y. <- ifelse(pies$roost_site == "Makira", pies$y. - 1, pies$y.)


pies$x. <- ifelse(pies$roost_site == "Ankarana", pies$x. + 2.5, pies$x.)
pies$y. <- ifelse(pies$roost_site == "Ankarana", pies$y. + 0.07, pies$y.)


head(pies)

#plot pie chart 
loko<-c("Rousettus madagascariensis"="#B200ED","Eidolon dupreanum"="#7FFF00","Pteropus rufus"="#0000FF")

p4 <- p2b+geom_path(data = pies, mapping = aes(x = x, y = y, group = roost_site), size = 0.25) + #tsy misy haiko aloha ny tena ilaina azy
  scale_size_identity() +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.8,.85),
                    legend.margin = margin(),
                    legend.title=element_text(face = 'bold', size = 8),
                    legend.background = element_rect(color="light grey",size = .1),
                    legend.text = element_text(size = 6,face = "italic"))+
  annotate("segment", x=pies$x, xend=pies$x.,y=pies$y,yend=pies$y.,size=.7,alpha=.5)+ # put the lines
  new_scale_fill() +
  geom_scatterpie(aes(x=x., y=y., group = roost_site, r = log10(nombre/20)), data = pies,
                  cols=c("bat_species"), long_format=TRUE, color = NA)+ 
  scale_fill_manual(values = loko,name="Bat species")+
  geom_scatterpie_legend(log10(c(44.959,50,200,1500)/25),
                         x=51.8, y=-23.5, 
                         n=3,
                         labeller = function(x) paste(((x)^2)*500,"indiv"))


#print(p4)

ggsave(file = paste0(homewd, "final-figures/Fig1B_final_map.pdf"),
       plot=p4,
       units="mm",  
       width=40, 
       height=60, 
       scale=3, 
       dpi=300)

