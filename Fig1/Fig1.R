rm(list=ls())

#packages
library(sf)
library(mapplots)
library(scatterpie)
library(maptools)
library(plyr) 
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
homewd = "/Users/carabrook/Developer/Mada-Bat-CoV/" 
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
  theme(plot.margin = unit(c(-1,.5,-1.5,.1),"cm"))+
  xlab("Longitude") + ylab("Latitude") 
#print(p1)
# # 
#   ggsave(file = paste0(homewd, "final-figures/tmp1.pdf"),
#          plot = p1,
#          units="mm",  
#          width=60, 
#          height=55, 
#          scale=3, 
#          dpi=300)
# 
  
  

#import CoV data
dat <- read.csv(file = paste0(homewd,"/metadata/all_NGS_8_3_2021_distribute.csv"), header = T, stringsAsFactors = F )
head(dat)
names(dat)

#only plot feces
dat = subset(dat, sample_type=="feces")

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
coordinate$label <- coordinate$species
coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum"
#load GPS point and label
p2b<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "#1B262C", size=3,
            nudge_x = c(-2,-3.6,8),
            nudge_y = c(3,-1.1,-.3),
            check_overlap = T)+
  annotation_scale(location = "bl", width_hint = 0.05) +    #scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.03, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)+
  geom_text_repel(segment.colour="black")+
  theme_bw() +theme(panel.grid = element_blank(), 
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
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

pies <- ddply(dat, .(species, roost_site, latitude_s, longitude_e, age_class, plot_class), summarise, value=length(sampleid))



tot_sum = ddply(pies,.(species, age_class), summarise,N=sum(value))

pies <- merge(pies, tot_sum, by=c("species", "age_class"), all.x=T)

pies$plot_class <- factor(pies$plot_class, levels=c( "juvenile: CoV neg", "adult: CoV neg", "juvenile: CoV pos", "adult: CoV pos"))

#now split into two pies
piesJ = subset(pies, age_class=="J")
piesA = subset(pies, age_class=="A")


###Get the pie data in the right format###
colz = c('adult: CoV neg' ="dodgerblue4", 'adult: CoV pos' ="firebrick4", 'juvenile: CoV neg' ="dodgerblue", 'juvenile: CoV pos' ="firebrick1")


p3<-ggplot() + 
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)), 
                  data = piesA, cols="plot_class", long_format=TRUE) +
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)), 
                  data = piesJ, cols="plot_class", long_format=TRUE) +
  scale_fill_manual(values=colz)


# 
# # copie of latitude (x.) and longitude (y.)
 piesA$x2 <- piesA$longitude_e
 piesA$y2 <- piesA$latitude_s
# 
# #manually move the pie chart in case there is an overlap (change x and y)
# 
 piesA$x2[piesA$species== "Pteropus rufus"] <- piesA$longitude_e[piesA$species== "Pteropus rufus"] -2
 piesA$y2[piesA$species== "Pteropus rufus"] <- piesA$latitude_s[piesA$species== "Pteropus rufus"] + 1
 
 
 piesA$x2[piesA$species== "Eidolon dupreanum"] <- piesA$longitude_e[piesA$species== "Eidolon dupreanum"] - .5
 piesA$y2[piesA$species== "Eidolon dupreanum"] <- piesA$latitude_s[piesA$species== "Eidolon dupreanum"] - 3
 
 piesA$x2[piesA$species== "Rousettus madagascariensis"] <- piesA$longitude_e[piesA$species== "Rousettus madagascariensis"] + 4
 piesA$y2[piesA$species== "Rousettus madagascariensis"] <- piesA$latitude_s[piesA$species== "Rousettus madagascariensis"] - 0
 
 head(piesA)
 
 
 # # copie of latitude (x.) and longitude (y.)
 piesJ$x2 <- piesJ$longitude_e
 piesJ$y2 <- piesJ$latitude_s
 # 
 # #manually move the pie chart in case there is an overlap (change x and y)
 # 
 piesJ$x2[piesJ$species== "Pteropus rufus"] <- piesJ$longitude_e[piesJ$species== "Pteropus rufus"] + 1
 piesJ$y2[piesJ$species== "Pteropus rufus"] <- piesJ$latitude_s[piesJ$species== "Pteropus rufus"] + 3
 
 
 piesJ$x2[piesJ$species== "Eidolon dupreanum"] <- piesJ$longitude_e[piesJ$species== "Eidolon dupreanum"] - 4.5
 piesJ$y2[piesJ$species== "Eidolon dupreanum"] <- piesJ$latitude_s[piesJ$species== "Eidolon dupreanum"] - 3
 
 piesJ$x2[piesJ$species== "Rousettus madagascariensis"] <- piesJ$longitude_e[piesJ$species== "Rousettus madagascariensis"] + 3
 piesJ$y2[piesJ$species== "Rousettus madagascariensis"] <- piesJ$latitude_s[piesJ$species== "Rousettus madagascariensis"] - 3
 
 head(piesJ)

#plot pie chart 
#loko<-c("Rousettus madagascariensis"="#B200ED","Eidolon dupreanum"="#7FFF00","Pteropus rufus"="#0000FF")

#this is Fig1A
p4 <- p2b+
  annotate("segment", x=piesA$longitude_e, xend=piesA$x2,y=piesA$latitude_s,yend=piesA$y2,size=.7)+ # put the lines
  annotate("segment", x=piesJ$longitude_e, xend=piesJ$x2,y=piesJ$latitude_s,yend=piesJ$y2,size=.7)+ # put the lines
  geom_scatterpie(aes(x=x2, y=y2, r=(log10(N)/1.2)), 
                          data = piesA, cols="plot_class", long_format=TRUE) +
  geom_scatterpie(aes(x=x2, y=y2, r=(log10(N)/1.2)), 
                  data = piesJ, cols="plot_class", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.8,.8),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 7.5)) +
  scale_fill_manual(values=colz) +
  geom_scatterpie_legend(log10(c(10,100)/1.2),
                         x=54.5, y=-23.5, 
                         n=2,
                         labeller = function(x) paste(10^(x)*1.2,"indiv"))

#print(p4)

Fig1a <- p4


#and this is Fig 1B

#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")

#check sites
unique(dat$roost_site) #all sites are in the moramanga area and can be treated as one

#now there are some urine and some fecal samples
#one of the urine samples is the same individual as one of the fecal samples
#so go ahead and remove the fecal sample (negative) from consideration



#get the date of the first day of every week
dat$epiwk <- cut(dat$collection_date, "week")
dat$epiwk <- as.Date(as.character(dat$epiwk))

#change CoV to numeric
dat$CoV[dat$CoV=="N"] <- 0
dat$CoV[dat$CoV=="Y"] <- 1
dat$CoV <- as.numeric(dat$CoV)

names(dat)[names(dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date

dat.list <- dlply(dat, .(sampleid))

#for those with multiple samples, slim to one


out = c(unlist(lapply(dat.list, nrow)))
out[out>1] 


slim.down <- function(df){
  if(nrow(df)==1){
    return(df)
  }else if (nrow(df)>1){
    max_CoV = max(df$CoV)
    df = df[1,]
    df$CoV = max_CoV
    return(df) 
  }
  
}

dat.list<- lapply(dat.list, slim.down)

dat.new <- data.table::rbindlist(dat.list)
out = c(unlist(lapply(dat.list, nrow)))
out[out>1] #none

dat <- dat.new
#should be two urine samples
length(dat$sample_type[dat$sample_type=="urine"])#2

#okay to go.

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species, age_class,epiwk), summarise, N=length(CoV), pos=sum(CoV))

#get negatives and prevalence
dat.sum$neg= dat.sum$N-dat.sum$pos
dat.sum$prevalence <- dat.sum$pos/dat.sum$N

#and confidence intervals on the prevalence
CIs <- mapply(FUN=prop.test, x=as.list(dat.sum$pos), n=as.list(dat.sum$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

#and extract the upper and lower CIs
get.CIs <- function(df){
  lci = df$conf.int[1]
  uci = df$conf.int[2]
  out.df<- cbind.data.frame(lci=lci, uci=uci)
  return(out.df)
}

CIs <- lapply(CIs, get.CIs)

dat.sum$lci <- c(unlist(sapply(CIs, '[',1)))
dat.sum$uci <- c(unlist(sapply(CIs, '[',2)))

#simplify= the name of "host_genus_species" 
names(dat.sum)[names(dat.sum)=="host_genus_species"] <- "species"

#and plot
#here's a vector assigning colors to each species
colz = c("Eidolon dupreanum"="steelblue1", "Pteropus rufus" = "violetred", "Rousettus madagascariensis" = "seagreen" )
names(dat.sum)[names(dat.sum)=="age_class"] <- "age"
dat.sum$age[dat.sum$age=="A"] <- "adult"
dat.sum$age[dat.sum$age=="J"] <- "juvenile"

shapez = c("juvenile" = 17, "adult" = 16)

p1 <- ggplot(data=dat.sum) + #here is the dataset
  geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species, group=age), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk, y= prevalence, color=species,shape=age, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("CoV Prevalence") + #change the name of the y-axis
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  scale_shape_manual(values=shapez) +
  facet_grid(age~.) +
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_blank())  #more style features


p1

#most information seems to be captured in plot 1, so lets keep that

#or try another plot that includes the winter season for Moramanga instead
seas.dat1 = cbind.data.frame(x=as.Date(c("2018-06-01","2018-09-01")), ymin=c(-1,-1), ymax=c(2,2), label="dry season")
seas.dat2 = cbind.data.frame(x=as.Date(c("2018-02-01","2018-06-01")), ymin=c(-1,-1), ymax=c(2,2), label="late-stage juveniles")

#jitter dates manually
dat.sum$epiwk_jitter <- dat.sum$epiwk
dat.sum$epiwk_jitter[dat.sum$species=="Rousettus madagascariensis"] <- dat.sum$epiwk_jitter[dat.sum$species=="Rousettus madagascariensis"] + 1
dat.sum$epiwk_jitter[dat.sum$species=="Eidolon dupreanum"] <- dat.sum$epiwk_jitter[dat.sum$species=="Eidolon dupreanum"] + 2

text.dat = cbind.data.frame(date=c(as.Date("2018-04-10"), as.Date("2018-07-15")),y=c(.97,.97), label=c("late-stage juveniles", "dry season"), age=c("juvenile", "juvenile"))

Fig1b <-  ggplot(data=dat.sum) +
  geom_ribbon(data=seas.dat1, aes(x=x, ymin=ymin, ymax=ymax), fill="cornflowerblue", alpha=.3) +
  geom_ribbon(data=seas.dat2, aes(x=x, ymin=ymin, ymax=ymax), fill="lightgoldenrod1", alpha=.3) +
  geom_errorbar(aes(x=epiwk_jitter, ymin=lci, ymax=uci, color=species,  group=age), size=.2) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk_jitter, y= prevalence, color=species, shape=age, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("CoV Prevalence") + #change the name of the y-axis
  geom_text(data=text.dat, aes(x=date,y=y, label=label), size=3) +
  #geom_text(data=text.dat, aes(x=as.Date("2018-04-24"),y=.95, label="late-stage juveniles"), size=3) +
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  scale_fill_manual(values=colz) + #assign the colors manually using the vector above
  scale_shape_manual(values=shapez) +#
  facet_grid(age~.) +
  theme_bw() + #some style features
  coord_cartesian(ylim=c(0,1)) +
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic", size = 8),
        legend.position = c(.28,.92),
        strip.background = element_rect(fill="white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill="white"),
        #legend.title = element_blank(),
        legend.spacing.y = unit(.05, "cm"),
        plot.margin = unit(c(.8,.2,.7,.3),"cm"),
        axis.title.x = element_blank())  #more style features

Fig1b



Fig1all <- cowplot::plot_grid(Fig1a, Fig1b, nrow=1, ncol=2, labels = c("(A)", "(B)"), label_x = -.01, label_y = .99)


ggsave(file = paste0(homewd, "final-figures/Fig1.pdf"),
       plot=Fig1all,
       units="mm",  
       width=150, 
       height=65, 
       scale=3, 
       dpi=300)




