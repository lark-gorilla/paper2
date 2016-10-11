# 10/10/16 Lancaster, UK
# Model HMM data modelling foraging against non-foraging


rm(list=ls())
library(ggplot2)
library(reshape2)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/GPS_LT_141516_hmm_oceano_attribs.csv", h=T, strip.white=T)

d1<-melt(dat[,c(1:5, 26:36)], id.vars=c("Latitude", "Longitude", "Date", "Time", "TrackID"))

p<-ggplot(d1, aes(variable, value))
p+geom_boxplot()+facet_wrap(~variable, scale="free")

# investigate and remove outliers

dat[which(dat$AsstAG>20),] # !!!! WRAP in which() otherwise you get loads of NA rows!!!
#for moment we'll just turn to NA
dat[which(dat$AsstAG>20),]$AsstAG<-NA 

dat[which(dat$chlMH1>0.3),]
# these are reef-inflated CHL values, hould probably leave rather than NA
# be aware in the heron model

qplot(chlMH1, data=dat[which(dat$chlMH1<0.3),], geom="histogram", bins=60)

#transformations
dat$LOGchlMH1<-log(dat$chlMH1)
dat$LOGchlVH3<-log(dat$chlVH3)
dat$SQRTsmt<-sqrt(dat$smt)

#collinearity
mcor<-cor(na.omit(dat[26:39]))
mcor 

library(car)
vif(lm(1:nrow(dat)~ekmU+modW+sstAG +sstOi+AsstAG+AsstOi+
         sshHy+bty+SQRTsmt+LOGchlVH3+LOGchlMH1,data=dat))

