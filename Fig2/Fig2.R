rm(list=ls())

library(ggplot2)   
library(plyr)
library(dplyr)
library(lubridate)
library(ISOweek)
library(cowplot)
library(mgcv)


#set wd()
#change for your computer
homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig2"))

#load fecal data
dat <- read.csv(file = paste0(homewd, "/metadata/bat_fecal_metadata_8_3_2021.csv"), header = T, stringsAsFactors = F)
head(dat)

#get into date form
dat$collection_date <- as.Date(dat$collection_date, format = "%m/%d/%y")

#check sites
unique(dat$site) #all sites are in the moramanga area and can be treated as one

#get the date of the first day of every week
dat$epiwk <- cut(dat$collection_date, "week")
dat$epiwk <- as.Date(as.character(dat$epiwk))

#change CoV to numeric
dat$CoV[dat$CoV=="N"] <- 0
dat$CoV[dat$CoV=="Y"] <- 1
dat$CoV <- as.numeric(dat$CoV)

#summarize into prevance by species and epiwk
dat.sum <- ddply(dat, .(host_genus_species, epiwk), summarise, N=length(CoV), pos=sum(CoV))

#get negatives and prevalence
dat.sum$neg= dat.sum$N-dat.sum$pos
dat.sum$prevalence <- dat.sum$pos/dat.sum$N

#and confidence intervals on the prevalence
CIs <- mapply(FUN=prop.test, x=as.list(dat.sum$pos), n=as.list(dat.sum$N), MoreArgs = list(alternative = "two.sided", conf.level = .95), SIMPLIFY = F)

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

p1 <- ggplot(data=dat.sum) + #here is the dataset
      geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
      geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
      ylab("CoV Prevalence") + #change the name of the y-axis
      scale_color_manual(values=colz) + #assign the colors manually using the vector above
      theme_bw() + #some style features
        theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
              axis.title.x = element_blank())  #more style features
      
      
p1


#also try plotting just by week of the year, independent of year

#get week of year
dat.sum$week <- isoweek(dat.sum$epiwk)

#summarise by week
dat.sum2<- ddply(dat.sum, .(species, week), summarise, N=sum(N), pos=sum(pos), neg=sum(neg))

dat.sum2$prevalence <- dat.sum2$pos/dat.sum2$N

CIs <- mapply(FUN=prop.test, x=as.list(dat.sum2$pos), n=as.list(dat.sum2$N), MoreArgs = list(alternative = "two.sided", conf.level = .95), SIMPLIFY = F)

CIs <- lapply(CIs, get.CIs)

dat.sum2$lci <- c(unlist(sapply(CIs, '[',1)))
dat.sum2$uci <- c(unlist(sapply(CIs, '[',2)))

#plot 

p2 <- ggplot(data=dat.sum2) + #here is the dataset
  geom_errorbar(aes(x=week, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=week, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("CoV Prevalence") + #change the name of the y-axis
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        axis.title.x = element_blank()) #more style features
  

p2

#and try ploting by month of year


#get week of year
dat.sum$month <- month(dat.sum$epiwk)

#summarise by week
dat.sum3<- ddply(dat.sum, .(species, month), summarise, N=sum(N), pos=sum(pos), neg=sum(neg))

dat.sum3$prevalence <- dat.sum3$pos/dat.sum3$N

CIs <- mapply(FUN=prop.test, x=as.list(dat.sum3$pos), n=as.list(dat.sum3$N), MoreArgs = list(alternative = "two.sided", conf.level = .95), SIMPLIFY = F)

CIs <- lapply(CIs, get.CIs)

dat.sum3$lci <- c(unlist(sapply(CIs, '[',1)))
dat.sum3$uci <- c(unlist(sapply(CIs, '[',2)))


p3 <- ggplot(data=dat.sum3) + #here is the dataset
  geom_errorbar(aes(x=month, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=month, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("CoV Prevalence") + #change the name of the y-axis
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        axis.title.x = element_blank()) + #more style features
  scale_x_continuous(breaks=1:12, labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

p3


#plot all three
cowplot::plot_grid( p1, p2,p3, nrow = 3, ncol=1) 


#most information seems to be captured in plot 1, so lets keep that

#Fit a GAM model to the seasonal prevalence, by species
#tom make these binomial models, first get the negatives at each time point

#subset each species dataset 
Pter.dat = subset(dat.sum, species=="Pteropus rufus")
Eid.dat = subset(dat.sum, species=="Eidolon dupreanum")
Rou.dat = subset(dat.sum, species=="Rousettus madagascariensis")


gam1 <- gam(cbind(pos,neg) ~ s(as.numeric(epiwk),k=7, bs="tp"),
                               family=binomial, data=Pter.dat)
summary(gam1)
#look at it
plot(gam1)

gam2 <- gam(cbind(pos,neg) ~ s(as.numeric(epiwk),k=6, bs="tp"),
            family=binomial, data=Eid.dat)
summary(gam2)
#look at it
plot(gam2)



gam3 <- gam(cbind(pos,neg) ~ s(as.numeric(epiwk),k=5, bs="tp"),
            family=binomial, data=Rou.dat)
summary(gam3)
#look at it
plot(gam3)


#not sure there is any point, 
#but can plot the modeled prevalence if desired

dat.sum$predicted_prevalence <- NA
dat.sum$predicted_prevalence_lci <- NA
dat.sum$predicted_prevalence_uci <- NA

dat.sum$predicted_prevalence[dat.sum$species=="Pteropus rufus"] <- predict.gam(gam1, type="response")
dat.sum$predicted_prevalence_lci[dat.sum$species=="Pteropus rufus"] <- dat.sum$predicted_prevalence[dat.sum$species=="Pteropus rufus"] - 1.96*predict.gam(gam1, type="response", se.fit = T)$se
dat.sum$predicted_prevalence_uci[dat.sum$species=="Pteropus rufus"] <- dat.sum$predicted_prevalence[dat.sum$species=="Pteropus rufus"] + 1.96*predict.gam(gam1, type="response", se.fit = T)$se

dat.sum$predicted_prevalence[dat.sum$species=="Eidolon dupreanum"] <- predict.gam(gam2, type="response")
dat.sum$predicted_prevalence_lci[dat.sum$species=="Eidolon dupreanum"] <- dat.sum$predicted_prevalence[dat.sum$species=="Eidolon dupreanum"] - 1.96*predict.gam(gam2, type="response", se.fit = T)$se
dat.sum$predicted_prevalence_uci[dat.sum$species=="Eidolon dupreanum"] <- dat.sum$predicted_prevalence[dat.sum$species=="Eidolon dupreanum"] + 1.96*predict.gam(gam2, type="response", se.fit = T)$se


dat.sum$predicted_prevalence[dat.sum$species=="Rousettus madagascariensis"] <- predict.gam(gam3, type="response")
dat.sum$predicted_prevalence_lci[dat.sum$species=="Rousettus madagascariensis"] <- dat.sum$predicted_prevalence[dat.sum$species=="Rousettus madagascariensis"] - 1.96*predict.gam(gam3, type="response", se.fit = T)$se
dat.sum$predicted_prevalence_uci[dat.sum$species=="Rousettus madagascariensis"] <- dat.sum$predicted_prevalence[dat.sum$species=="Rousettus madagascariensis"] + 1.96*predict.gam(gam3, type="response", se.fit = T)$se

#and plot with the data

p4 <-  ggplot(data=dat.sum) +
        geom_ribbon(aes(x=epiwk, ymin=predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=species), alpha=.2) +
        #geom_line(data=dat.sum, aes(x=epiwk, y=predicted_prevalence, color=species), alpha=.3) + 
        geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
        geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
        ylab("CoV Prevalence") + #change the name of the y-axis
        scale_color_manual(values=colz) + #assign the colors manually using the vector above
        scale_fill_manual(values=colz) + #assign the colors manually using the vector above
        theme_bw() + #some style features
        theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        axis.title.x = element_blank()) #more style features
p4


#or try another plot that includes the winter season for Moramanga instead
seas.dat = cbind.data.frame(x=as.Date(c("2018-06-01","2018-09-01")), ymin=c(0,0), ymax=c(1,1))

p5 <-  ggplot(data=dat.sum) +
  geom_ribbon(data=seas.dat, aes(x=x, ymin=ymin, ymax=ymax), fill="cornflowerblue", alpha=.3) +
  geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("CoV Prevalence") + #change the name of the y-axis
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  scale_fill_manual(values=colz) + #assign the colors manually using the vector above
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        axis.title.x = element_blank()) #more style features
p5

#or combine the two


p6 <-  ggplot(data=dat.sum) +
       geom_vline(data=seas.dat, aes(xintercept=x), color="navy", size=1, linetype=2) +
       geom_ribbon(aes(x=epiwk, ymin=predicted_prevalence_lci, ymax=predicted_prevalence_uci, fill=species), alpha=.2) +
       #geom_line(data=dat.sum, aes(x=epiwk, y=predicted_prevalence, color=species), alpha=.3) + 
       geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
       geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
       ylab("CoV Prevalence") + #change the name of the y-axis
       scale_color_manual(values=colz) + #assign the colors manually using the vector above
       scale_fill_manual(values=colz) + #assign the colors manually using the vector above
       theme_bw() + #some style features
       coord_cartesian(ylim=c(0,1)) + #set ylimits
       theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
      axis.title.x = element_blank()) #more style features
p6

#and save plot