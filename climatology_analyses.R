# 25/09/16 Lancaster, UK
# Model climatology data comparing foraging areas against background 

rm(list=ls())
library(ggplot2)
library(reshape2)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/paper2_extractionV2.csv", h=T)

d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "dtyp"))

p<-ggplot(d1, aes(dset, value))
p+geom_boxplot()+facet_wrap(~variable, scale="free")

## We're going to use the Information Theoretic approach to test presence (50%UD) absence (max range)
## of WTSH using a ?glm? testing oceanographic vs tuna variables using AIC. We will do this 
## seperately for each colony.

# first we need to remove outliers and transform some variables

qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")

qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

# simplist way is do an na.omit, which basically kills all coastline pixels
qplot(Longitude, Latitude, data=dat)+geom_point(data=na.omit(dat), aes(Longitude, Latitude), colour=2)

dat<-na.omit(dat)
#now how does it look?
d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "dtyp"))
qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")
qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

# seamounts and chl still need some attention
qplot(Longitude, Latitude, data=dat)+geom_point(data=dat[dat$chl<1000,], aes(Longitude, Latitude), colour=2)
## hmm although these are likely reef effects its better to leave in and transform
qplot(chl, data=dat);qplot(log(chl), data=dat)



