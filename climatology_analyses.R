# 25/09/16 Lancaster, UK
# Model climatology data comparing foraging areas against background 

rm(list=ls())
library(ggplot2)
library(reshape2)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/paper2_extractionV2.csv", h=T, strip.white=T)

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
# bit too brutal.
# best plan is too set extreme chl values to NA then use those to mask out other suspect variables

qplot(Longitude, Latitude, data=dat)+geom_point(data=dat[dat$chl<1500,], aes(Longitude, Latitude), colour=2)

dat[!is.na(dat$chl) & dat$chl>1500,]$chl<-NA
# also knock out some dodgy looking wind values
dat[is.na(dat$chl),]$wnd<-NA
# and some weird sst
dat[!is.na(dat$sst) & dat$sst<0,]$sst<-NA


# how we looking now
d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "dtyp"))
qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")
qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

g<-ggplot(data=d1, aes(x=factor(variable), y=value))
g+geom_violin(alpha=0.5, color="gray")+
  geom_jitter(alpha=0.5,position = position_jitter(width = 0.1))+
  facet_wrap(~variable, scale="free")

#weird ekem point
dat[!is.na(dat$ekm) &dat$ekm< -0.00001,]$ekm<-NA

# seamounts and chl still need some attention
qplot(chl, data=dat);qplot(log(chl), data=dat)
qplot(smt, data=dat);qplot(sqrt(smt), data=dat)

dat$chl_log<-log(dat$chl)
dat$smt_sqt<-sqrt(dat$smt)
dat$PA<-0
dat[dat$dtyp=="ud50_pres",]$PA<-1 # making presence absence variable
dat_heron<-dat[dat$dset=="heronBuff"  | dat$dset=="LTHeronGPS2015"  | 
                 dat$dset=="LTHeronPTT2011"  | dat$dset=="LTHeronPTT2013" , ]

dat_lhi<-dat[dat$dset=="lhiBuff"  | dat$dset=="LTLHIGPS2014"  | 
                 dat$dset=="LTLHIGPS2015"  | dat$dset=="LTLHIGPS2016" , ]

# Aight lets do some modellin'

# check for multicollinearoty
mcor<-cor(na.omit(enviro[,-c(1:4)]))
mcor # seems to be some collinearity
library(car)
vif(lm(1:nrow(dat_heron)~sst+shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+bet_tot_03_ave+yft_adu_03_ave+
         yft_juv_03_ave+yft_tot_03_ave+skj_adu_03_ave+
         skj_juv_03_ave+skj_tot_03_ave, data=dat_heron))

vif(lm(1:nrow(dat_lhi)~sst+shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+bet_tot_03_ave+yft_adu_03_ave+
         yft_juv_03_ave+yft_tot_03_ave+skj_adu_03_ave+
         skj_juv_03_ave+skj_tot_03_ave, data=dat_lhi))

# dropping sst, yft_tot, bet_tot, skj_tot, yft_juv and skj_juv seems to drop all below VIF 10

vif(lm(1:nrow(dat_lhi)~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+yft_adu_03_ave+
         skj_adu_03_ave
         , data=dat_lhi))

vif(lm(1:nrow(dat_heron)~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+yft_adu_03_ave+
         skj_adu_03_ave
       , data=dat_heron))

# ok looks good, all values < 10

#library(PerformanceAnalytics)
#chart.Correlation(dat_heron[,c(6,7,9,10,12,13,14,16,18,22)], histogram=F, pch=19)

mcor<-cor(na.omit(dat_heron[,-c(1:4)]))
mcor # we can see SKJ adu (0.8)and YFT adu (0.67) are highly coorelated with their juvs
mcor<-cor(na.omit(dat_lhi[,-c(1:4)]))
mcor #we can see SKJ adu (0.8) is highly coorelated with its juv, YFT less so (0.2)

### try some models

m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          bet_juv_03_ave+yft_adu_03_ave+skj_adu_03_ave,
          data=dat_heron, family="binomial")
summary(m1)
sum((resid(m1, type="pearson")^2))/df.residual(m1)

library(MuMIn)
options(na.action = "na.fail") 
dr<-dredge(m1)

# If i run each year seperately the model becomes zero inflated. I can either sub
# sample the absence data to reduce the zinf then run individually or run with
# a random effect for dset. I guess i need to decide if i want to ask if the pattern
# varies between year, ie is it consistent or just the general pattern for each col... maybe the latter

