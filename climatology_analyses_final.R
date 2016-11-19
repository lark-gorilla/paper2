# 19/11/16 LBishops Castle, UK
# Model climatology data comparing foraging areas against background 
# Final paper version including method to deal with SPAC

rm(list=ls())
library(ggplot2)
library(reshape2)

# extra functions for pairs multicollinearity exploration
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

## read in data

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/paper2_extractionV5.csv", h=T, strip.white=T)
# V3 takes all PA points

d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "Count"))

p<-ggplot(d1, aes(dset, value))
p+geom_boxplot()+facet_wrap(~variable, scale="free")

## We're going to use model selection to test presence (50%UD) absence (max range)
## of WTSH using a glm testing oceanographic vs tuna variables using AIC. We will do this 
## seperately for each colony.

# first we need to remove outliers and transform some variables

qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")

qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

qplot(Longitude, Latitude, data=dat)+geom_point(data=dat[dat$chl<1500,], aes(Longitude, Latitude), colour=2)

dat[which(dat$chl>1500),]$chl<-NA
# also knock out some dodgy looking wind values
dat[is.na(dat$chl),]$wnd<-NA
# and some weird sst
dat[which(dat$sst<0),]$sst<-NA
#weird ekem point
dat[which(dat$ekm< -0.00001),]$ekm<-NA
#remove skj_juv outliers
dat[which(dat$skj_juv_03_ave>28),]$skj_juv_03_ave<-NA


# seamounts and chl still need some attention
qplot(chl, data=dat);qplot(log(chl), data=dat)
qplot(smt, data=dat);qplot(sqrt(smt), data=dat)

dat$chl_log<-log(dat$chl)
dat$smt_sqt<-sqrt(dat$smt)
dat$PA<-0
dat[dat$Count>0,]$PA<-1 # making presence absence variable

# remove some outlier presence clusters that cause modeeling issues later
dat<-dat[-which(dat$dset=="Heron"& dat$PA==1 & dat$Longitude<153.04 &
                  dat$Latitude< -21.44),]

dat<-dat[-which(dat$dset=="Heron"& dat$PA==1 & dat$Longitude<152.5 &
                  dat$Latitude< -18.74),]

dat<-dat[-which(dat$dset=="Heron" & dat$PA==1 & dat$Longitude<151.4),]

## NOW we do naomit! 
## This is integral to modelling with the IT approach and others, because
## when subset models are created, they have different numbers of NAs
## and therefore data: they are not comparable

nrow(dat)
nrow(na.omit(dat))
dat<-na.omit(dat)


# how we looking now
d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "Count", "PA"))
qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")
qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

#g<-ggplot(data=d1, aes(x=factor(variable), y=value))
#g+geom_violin(alpha=0.5, color="gray")+
#  geom_jitter(alpha=0.5,position = position_jitter(width = 0.1))+
#  facet_wrap(~variable, scale="free")

# setup colony datasets


####### Heron Island analyses #######

heronP<-dat[dat$dset=="Heron" & dat$Count>0,]
heronA<-dat[dat$dset=="Heron" & dat$Count==0,]
heronPA<-heronA[sample(1:nrow(heronA), 3200),]

dat_heron<-rbind(heronP, heronPA)

# Below code removed duplicate points where both presence and absence occur
# together, however caused perfect seperation in binomial glm
#nrow(dat_heron)   
#sp1<-SpatialPointsDataFrame(dat_heron[,1:2], data=dat_heron)
#sp2<-remove.duplicates(sp1, remove.second = TRUE, memcmp = TRUE) # removes duplicate points, remove.second makes sure presences remain and absences are removed
#dat_heron<-sp2@data
#nrow(dat_heron) 

lhiP<-dat[dat$dset=="LHI" & dat$Count>0,]
lhiA<-dat[dat$dset=="LHI" & dat$Count==0,]
lhiPA<-lhiA[sample(1:nrow(lhiA), 3200),]
dat_lhi<-rbind(lhiP, lhiPA)

# ok attempt with ME
library(spdep)

sp1<-SpatialPointsDataFrame(dat_heron[,1:2], data=dat_heron)

dat.nb4<-knearneigh(sp1@coords, k=3, longlat=T)

dat.nb4<-knn2nb(dat.nb4)
dat.nb4<-make.sym.nb(dat.nb4)
dat.wt4<-nb2listw(dat.nb4, style="W")

plot(sp1, col = "grey60")
plot(dat.nb4, coordinates(sp1), pch = 19, cex = 0.3, add = TRUE)

m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+
          bet_juv_03_ave+yft_adu_03_ave+
          skj_juv_03_ave, data=dat_heron, 
        family=binomial)

me.fit<-ME(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+
             bet_juv_03_ave+yft_adu_03_ave+
             skj_juv_03_ave, data=dat_heron, 
           family=binomial, listw = dat.wt4, verbose=T,alpha=0.05 )


# z-transform (scale and center) to make variables comparable on same scale
# not doing that now!
#library(vegan)
#dat_heron<-cbind(dat_heron[,c(1,2,24:25)],
#                 decostand(dat_heron[,c(5:7, 9,10,12:23)], method="standardize"))


# remove some other problematic datapoints not now using interpolated tuna data
#dat_heron<-dat_heron[-which(row.names(dat_heron)%in%c("72622", "72614", "72606")),]

# now have a look at collinearity
library(car)
vif(lm(1:nrow(dat_heron)~sst+shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+bet_tot_03_ave+yft_adu_03_ave+
         yft_juv_03_ave+yft_tot_03_ave+skj_adu_03_ave+
         skj_juv_03_ave+skj_tot_03_ave, data=dat_heron))  

vif(lm(1:nrow(dat_heron)~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+yft_adu_03_ave+
         skj_adu_03_ave,
       data=dat_heron))

# ok looks good, all values < 10, hmm skjadu and yftadu are corr according to 
cor(dat_heron[,-(1:4)])

#pairs(dat_heron[,-(1:2)], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_heron[,11:19], upper.panel = panel.smooth,lower.panel=panel.cor)
# ok so actually for Heron I include both BET then
# adult yft and juv skj
# final one for model
pairs(dat_heron[,c(6:12,15,17, 20, 21) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# looking at it in full, i would be tempted to remove shd and yellowfin adult


m1<-glm(Count~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          bet_juv_03_ave+yft_adu_03_ave+skj_juv_03_ave,
        data=dat_heron, family="poisson")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_heron$PA, fitted(m1))$A)

# now we see which variables are better suited to a second deg polynomial

d1<-melt(dat_heron, id.vars=c("PA", "Count", "dset"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=PA, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")

# we would do below for each variable individually to see the suitability of a poly

m_lin<-glm(PA~yft_adu_03_ave,data=dat_heron, family="binomial")
m_pol<-glm(PA~poly(yft_adu_03_ave, 2),data=dat_heron, family="binomial")
m_gam<-gam(PA~s(yft_adu_03_ave, bs="cr", k=3),data=dat_heron, family="binomial")
d1<-data.frame(PA=dat_heron$PA, yft_adu_03_ave=dat_heron$yft_adu_03_ave, m_lin=fitted(m_lin), m_pol=fitted(m_pol), m_gam=fitted(m_gam))
g1<-ggplot(data=d1, aes(y=PA, x=yft_adu_03_ave))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=m_lin, x=yft_adu_03_ave), colour=2)+
  geom_line(aes(y=m_pol, x=yft_adu_03_ave), colour=3)+geom_line(aes(y=m_gam, x=yft_adu_03_ave), colour=4)
anova(m_lin, m_pol, m_gam)
#poly model better in the case of bet_adu it gives the 0 to 1 warning so only use linear

# It was important to work out what was causing the warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
# it was the  BET_ADU poly so its now removed and we can proceed!

### RAC

m2<-glm(PA~poly(shd,2)+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          bet_juv_03_ave+poly(skj_juv_03_ave,2),
        data=dat_heron, family="binomial") 

m2<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+
          bet_adu_03_ave+yft_adu_03_ave+
          skj_juv_03_ave,
        data=dat_heron, family="binomial") 

resglm<-residuals(m2, type="pearson")

library(ncf)
corglm <- correlog(dat_heron$Longitude, dat_heron$Latitude, residuals(m2, type="pearson"), na.rm=T,
                   latlon=T, increment=30,resamp=1)

library(spdep)
RAC<-autocov_dist(resglm, cbind(dat_heron[,1], dat_heron[,2]),
                  nbs = 50, type = "inverse", zero.policy = T,
                  style = "B", longlat=TRUE)
# look paratem

makeRAC<-function(mod=glmX, ras=rasTempl, dat=trfdf){  
  values(ras)<-NA
  xy_residuals <-cbind(dat$Longitude, dat$Latitude, residuals(mod, type="pearson"))
  ras[cellFromXY(ras,xy_residuals[,1:2])]<-xy_residuals[,3]
  focal_rac_rast<-focal(ras, w=matrix(1,3,3), fun = mean, na.rm = TRUE)
  focal_rac_vect<-extract(focal_rac_rast,xy_residuals[,1:2])
  return(focal_rac_vect)}

library(raster)
chl<-raster("~/grive/phd/sourced_data/env_data/climatologies/CPMbfp12016-03-01.nc")

RAC2<-makeRAC(m2, chl, dat_heron)

dat_heron<-cbind(dat_heron, RAC2)
dat_heron<-cbind(dat_heron, RAC)

m3<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+
          bet_adu_03_ave+yft_adu_03_ave+
          skj_juv_03_ave+RAC,
        data=dat_heron, family="binomial") # nope

corglm2 <- correlog(dat_heron$Longitude, dat_heron$Latitude, residuals(m3, type="pearson"), na.rm=T,
                    latlon=T, increment=20,resamp=1)

library(arm)
m4<-bayesglm(PA~poly(shd,2)+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
               bet_adu_03_ave+bet_juv_03_ave+
               poly(skj_juv_03_ave,2)+RAC,
             data=dat_heron, family="binomial")


