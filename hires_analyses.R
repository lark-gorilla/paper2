# 10/10/16 Lancaster, UK
# Model HMM data modelling foraging against non-foraging


rm(list=ls())
library(ggplot2)
library(reshape2)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/GPS_141516_LT_hmm_oceano_attribs.csv", h=T, strip.white=T)

d1<-melt(dat[,c(1:5, 26:34)], id.vars=c("Latitude", "Longitude", "Date", "Time", "TrackID"))

p<-ggplot(d1, aes(value))
p+geom_boxplot()+facet_wrap(~variable, scale="free")
