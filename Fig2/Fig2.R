rm(list=ls())

library(ggplot2)   
library(plyr)
library(dplyr)
library(lubridate)
library(ISOweek)
library(cowplot)
library(mgcv)
library(lmodel2)


#load data 
homewd= "/Users/caraebrook/Documents/R/R_repositories/Mada-Bat-CoV"
setwd(paste0(homewd, "/Fig2"))
dat <- read.csv(file=paste0(homewd, "/metadata/all_NGS_8_3_2021_distribute.csv"), header = T, stringsAsFactors = F)


#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")

#check sites
unique(dat$roost_site) #all sites are in the moramanga area and can be treated as one

#now there are some urine and some fecal samples

#get the date of the first day of every week
dat$epiwk <- cut(dat$collection_date, "week")
dat$epiwk <- as.Date(as.character(dat$epiwk))

#change CoV to numeric
dat$CoV[dat$CoV=="N"] <- 0
dat$CoV[dat$CoV=="Y"] <- 1
dat$CoV <- as.numeric(dat$CoV)

names(dat)[names(dat)=="bat_species"] <- "species"

length(dat$CoV[dat$CoV==1 & dat$sample_type=="urine"])
#Only 2 urine positives. We assume these represent contamination of feces in urine.

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
dat.sum <- ddply(dat, .(species, epiwk), summarise, N=length(CoV), pos=sum(CoV))

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

CIs <- mapply(FUN=prop.test, x=as.list(dat.sum2$pos), n=as.list(dat.sum2$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

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


#get week of year
dat.sum$month <- month(dat.sum$epiwk)

#summarise by week
dat.sum3<- ddply(dat.sum, .(species, month), summarise, N=sum(N), pos=sum(pos), neg=sum(neg))

dat.sum3$prevalence <- dat.sum3$pos/dat.sum3$N

CIs <- mapply(FUN=prop.test, x=as.list(dat.sum3$pos), n=as.list(dat.sum3$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

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
summary(gam2) #rubbish
#look at it
plot(gam2)



gam3 <- gam(cbind(pos,neg) ~ s(as.numeric(epiwk),k=5, bs="tp"),
            family=binomial, data=Rou.dat)
summary(gam3)
#look at it
plot(gam3)


#not sure there is any point, because the model fits are so bad
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
seas.dat = cbind.data.frame(x=as.Date(c("2018-06-01","2018-09-01")), ymin=c(-1,-1), ymax=c(2,2), label="dry season")

p5 <-  ggplot(data=dat.sum) +
  geom_ribbon(data=seas.dat, aes(x=x, ymin=ymin, ymax=ymax), fill="cornflowerblue", alpha=.3) +
  geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species, group=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("CoV Prevalence") + #change the name of the y-axis
  geom_text(aes(x=as.Date("2018-07-15"),y=.95, label="dry season"), size=3) +
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  scale_fill_manual(values=colz) + #assign the colors manually using the vector above
  theme_bw() + #some style features
  coord_cartesian(ylim=c(0,1)) +
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        axis.title.x = element_blank()) #more style features
p5

#or combine the two


p6 <-  ggplot(data=dat.sum) +
       geom_ribbon(data=seas.dat, aes(x=x, ymin=ymin, ymax=ymax), 
                   color="black", fill="gray",alpha=.2) +
       geom_text(aes(x=as.Date("2018-07-15"),y=.95, label="dry season"), size=3) +
       #geom_vline(data=seas.dat, aes(xintercept=x), color="black", size=.7, linetype=2) +
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



#plot prevalence by age class and sex

#clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat$young_of_year)
dat$age_class <- dat$bat_age_class
dat$age_class[dat$age_class=="P" | dat$age_class=="L"] <- "A"
dat$age_class[dat$age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$age_class[dat$young_of_year=="yes"] <- "J"

#by age class
prev.dat <- ddply(dat, .(species, bat_sex, age_class), summarise, pos=sum(CoV), N=length(CoV))
prev.dat$prevalence = prev.dat$pos/prev.dat$N

prev.dat$age_class[prev.dat$age_class=="J"] <- "juveniles"
prev.dat$age_class[prev.dat$age_class=="A"] <- "adults"

prev.dat$age_class <- factor(prev.dat$age_class, levels=c("juveniles", "adults"))


#and confidence intervals on the prevalence
CIs <- mapply(FUN=prop.test, x=as.list(prev.dat$pos), n=as.list(prev.dat$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

#and extract the upper and lower CIs
get.CIs <- function(df){
  lci = df$conf.int[1]
  uci = df$conf.int[2]
  out.df<- cbind.data.frame(lci=lci, uci=uci)
  return(out.df)
}

CIs <- lapply(CIs, get.CIs)

prev.dat$lci <- c(unlist(sapply(CIs, '[',1)))
prev.dat$uci <- c(unlist(sapply(CIs, '[',2)))

#seems to be higher prevalence in juveniles -- could be driving transmission
colz = c("female" ="pink", "male" = "cornflowerblue")
p7 <- ggplot(data=prev.dat, aes(x=age_class, y=prevalence, fill=bat_sex)) + 
  scale_fill_manual(values=colz) +
  geom_bar( stat = "identity", position =position_dodge(.9) ) +
  facet_grid(~species) + theme_bw() + 
  geom_errorbar(aes(ymin=lci, ymax=uci), color="gray50", size=.5, width=.1, position = position_dodge(.9)) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        strip.text = element_text(face="italic"), legend.position = c(.9,.8), 
        legend.title = element_blank(), axis.title.x = element_blank())
p7


ggsave(file = paste0(homewd, "/final-figures/Fig2_seasonal_CoV_prev_by_species.png"),
       plot = p6,
       units="mm",  
       width=70, 
       height=40, 
       scale=3, 
       dpi=300)



########################################################
########################################################

#also check any recaps
head(dat)
dat$unique_tag_id[dat$unique_tag_id==""] <- NA
dat$unique_tag_id[is.na(dat$unique_tag_id)] <- dat$sampleid[is.na(dat$unique_tag_id)] 
dat.recap.sum <- ddply(dat, .(unique_tag_id), summarise, N=length(CoV), recap=unique(recap))
dat.recap.sum = subset(dat.recap.sum, N>1) #these are those that have been sampled more than once
#8 individuals

#edit over in main dataset
dat$recap_CoV <- "no"
for(i in 1:length(dat.recap.sum$unique_tag_id)){
  dat$recap_CoV[dat$unique_tag_id==dat.recap.sum$unique_tag_id[i]]  <- "yes"
}
dat.recap = subset(dat, recap_CoV=="yes") #16. 8 individuals
length(unique(dat.recap$unique_tag_id)) #8

#now plot as a timeline
dat.recap <- arrange(dat.recap, unique_tag_id, collection_date)
head(dat.recap)

#split and supply a different number for each
recap.list <- dlply(dat.recap, .(unique_tag_id))

recap.rearrange <- function(df){
  
    print(df$unique_tag_id)
    start=min(df$collection_date)
    end = max(df$collection_date)
    CoV_start = df$CoV[df$collection_date==start] 
    CoV_end = df$CoV[df$collection_date==end] 
    
    df <- dplyr::select(df, species, bat_sex, unique_tag_id)
    df <- df[!duplicated(df),]
    df$x = start
    df$CoV_start = CoV_start
    df$xend = end
    df$CoV_end= CoV_end
    return(df)
 
}

out <- lapply(recap.list, recap.rearrange)

plot.recap <- data.table::rbindlist(out)
head(plot.recap)

#plot each as a line segment
colz = c("Eidolon dupreanum"="steelblue1", "Pteropus rufus" = "violetred", "Rousettus madagascariensis" = "seagreen" )



p10 <- ggplot(data=plot.recap) + 
       geom_point(aes(x=x, y=CoV_start, color=species), shape=16) +
       geom_point(aes(x=xend, y=CoV_end, color=species), shape=17) +
       geom_segment(aes(x=x, y=CoV_start, xend=xend, yend=CoV_end, color=species)) +
       scale_color_manual(values=colz) + theme_bw() +
       theme(panel.grid = element_blank(), axis.title.x = element_blank())
  
p10

#most bats just stayed CoV negative, but two Eidolon progressed from negative to positive
#neg= Jul 2018, pos = Jan 2019
#neg= Feb 2018, pos = Apr 2018

#save the table
write.csv(plot.recap, file = "recap_table.csv", row.names = F)
