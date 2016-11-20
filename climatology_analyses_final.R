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

#write.csv(dat, "test4.csv", quote=F, row.names=F) # write out na omitted dat


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

#going for 1:3 presence-absence this geives models enough
# absence locations to capture background but not too many
# to supress the patterns in prediction plots (ie 10 zeros to every 1,
# means very little observable effect)

heronP<-dat[dat$dset=="Heron" & dat$Count>0,]
heronA<-dat[dat$dset=="Heron" & dat$Count==0,]
heronPA<-heronA[sample(1:nrow(heronA), 3000),]
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
lhiPA<-lhiA[sample(1:nrow(lhiA), 10000),]
dat_lhi<-rbind(lhiP, lhiPA)

# ok attempt with ME
library(spdep)

sp1<-SpatialPointsDataFrame(dat_heron[,1:2], data=dat_heron, proj4string=CRS("+proj=longlat + datum=wgs84"))

dat.nb4<-knearneigh(sp1@coords, k=3, longlat=T)
dat.nb4<-knn2nb(dat.nb4)
dat.nb4<-make.sym.nb(dat.nb4)
dat.wt4<-nb2listw(dat.nb4, style="W")
plot(sp1, col = "grey60")
plot(w, coordinates(sp1), pch = 19, cex = 0.3, add = TRUE)

#using distance

nb <- dnearneigh(sp1@coords, 10,30, longlat=T) 
w <- nb2listw(nb,style="B",zero.policy=T)
plot(sp1, col = "grey60")
plot(w, coordinates(sp1), pch = 19, cex = 0.3, add = TRUE)

me.fit<-ME(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+
             bet_juv_03_ave+yft_adu_03_ave+
             skj_juv_03_ave, data=dat_heron, 
           family=binomial, listw = dat.wt4, verbose=T,alpha=0.05 )


# finding distance to nearest point after na.omit, this will identify isolated points
library(sp)
library(rgeos)
# proj dat into m units for gDistance
sp1<-SpatialPointsDataFrame(dat_heron[,1:2], data=dat_heron, proj4string=CRS("+proj=longlat + datum=wgs84"))
sp1proj <- spTransform(sp1, CRS=CRS("+proj=laea +lon_0=155 +lat_0=-21"))

d <- gDistance(sp1proj, byid=T)

min.d <- apply(d, 1, function(x) min(x[x>1000])) # pulls nearest point, which isnt itself or within 1000m projection error range
dat_heron$min.d<-min.d
hist(dat_heron$min.d)

plot(Latitude~Longitude, dat_heron)
points(Latitude~Longitude, dat_heron[dat_heron$min.d>26000,], col=2)

dat_heron<-dat_heron[-which(dat_heron$min.d>26000),] #remove poor neighbors!

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
pairs(dat_heron[,c(6,7,9,10,12,13,14,17, 19, 22, 23) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# looking at it in full, i would be tempted to remove shd, tmc and yellowfin adult


m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          bet_juv_03_ave+yft_adu_03_ave+skj_juv_03_ave,
        data=dat_heron, family="binomial")
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

m_lin<-glm(PA~skj_juv_03_ave,data=dat_heron, family="binomial")
m_pol<-glm(PA~poly(skj_juv_03_ave, 2),data=dat_heron, family="binomial")
d1<-data.frame(PA=dat_heron$PA, skj_juv_03_ave=dat_heron$skj_juv_03_ave, m_lin=fitted(m_lin), m_pol=fitted(m_pol))
g1<-ggplot(data=d1, aes(y=PA, x=skj_juv_03_ave))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=m_lin, x=skj_juv_03_ave), colour=2)+
  geom_line(aes(y=m_pol, x=skj_juv_03_ave), colour=3)
anova(m_lin, m_pol)
#poly model better in the case of skj_juv it gives the 0 to 1 warning so only use linear

# It was important to work out what was causing the warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
# it was the  skj_juv poly so its now removed and we can proceed!

### RAC

m2<-glm(PA~shd+ekm+chl_log+wnd+poly(tmc,2)+
          smt_sqt+poly(bty,2)+
          bet_adu_03_ave+bet_juv_03_ave+
          poly(yft_adu_03_ave,2)+skj_juv_03_ave,
        data=dat_heron, family="binomial") 
summary(m2)
print(sum((resid(m2, type="pearson")^2))/df.residual(m2))
library(verification)
print(roc.area(dat_heron$PA, fitted(m2))$A)

resglm<-residuals(m2, type="pearson")
sp2<-SpatialPointsDataFrame(dat_heron[,1:2], data=data.frame(resglm), 
    proj4string=CRS("+proj=longlat + datum=wgs84"))
bubble(sp2, zcol='resglm')

#library(ncf)
#corglm <- correlog(dat_heron$Longitude, dat_heron$Latitude, residuals(m2, type="pearson"), na.rm=T,
#                   latlon=T, increment=25,resamp=1)
library(spdep)
RAC<-autocov_dist(resglm, cbind(dat_heron[,1], dat_heron[,2]),
                  nbs = 25, type = "one", zero.policy = T,
                  style = "B", longlat=TRUE)
# look paratem

dat_heron$RAC<-RAC

m3<-glm(PA~shd+ekm+chl_log+wnd+poly(tmc,2)+
          smt_sqt+poly(bty,2)+bet_adu_03_ave+bet_juv_03_ave+
          poly(yft_adu_03_ave,2)+skj_juv_03_ave+RAC,
        data=dat_heron, family="binomial")

summary(m3)
print(sum((resid(m3, type="pearson")^2))/df.residual(m3))
resglm<-residuals(m3, type="pearson")
sp2<-SpatialPointsDataFrame(dat_heron[,1:2], data=data.frame(resglm), 
                            proj4string=CRS("+proj=longlat + datum=wgs84"))
bubble(sp2, zcol='resglm')
plot(m3)

#plots of resids show some outliers messing with the program, find em and remove then refit

dat_heron<-dat_heron[-which(resglm< -25),]

m2<-glm(PA~shd+ekm+chl_log+wnd+poly(tmc,2)+
          smt_sqt+poly(bty,2)+
          bet_adu_03_ave+bet_juv_03_ave+
          poly(yft_adu_03_ave,2)+skj_juv_03_ave,
        data=dat_heron, family="binomial") 

RAC<-autocov_dist(residuals(m2, type="pearson"), cbind(dat_heron[,1], dat_heron[,2]),
                  nbs = 25, type = "one", zero.policy = T,
                  style = "B", longlat=TRUE)

dat_heron$RAC<-RAC

m3<-glm(PA~shd+ekm+chl_log+wnd+poly(tmc,2)+
          smt_sqt+poly(bty,2)+bet_adu_03_ave+bet_juv_03_ave+
          poly(yft_adu_03_ave,2)+skj_juv_03_ave+RAC,
        data=dat_heron, family="binomial")

plot(m3)
print(sum((resid(m3, type="pearson")^2))/df.residual(m3))
resglm<-residuals(m3, type="pearson")
print(roc.area(dat_heron$PA, fitted(m3))$A)
library(pscl)
pR2(m3)
summary(m3)

# looking at variable importance/contribution
summary(m3) # look to drop ekm and yrid and probs yft
anova(m3) # anova tests terms sequentially
library(survey)
regTermTest(m3, "ekm")
anova(m3, update(m3, ~.- ekm)) # basically tells the same as anova but gives significance too

library(car)
Anova(m3) # Anova from car tests terms according to marginality
# i.e. after all other terms have been included
# not much support from ekm, or skj

# See if step methods give same result from the m1 model

for.aic <- step(glm(PA~1,data=dat_heron, family="binomial"),
                direction = "forward", scope = formula(m3), k = 2, trace = 1) # forward AIC
for.bic <- step(glm(PA~1,data=dat_heron, family="binomial"),
                direction = "forward", scope = formula(m3), k = log(nrow(dat_heron)), trace = 0) # forward BIC
back.aic <- step(m3, direction = "backward", k = 2, trace = 0) # backward AIC
back.bic <- step(m3, direction = "backward", k = log(nrow(dat_heron)), trace = 0) # backward BIC

formula(for.aic);formula(for.bic);formula(back.aic);formula(back.bic);

formula(m2)
# so basically m2 and aic methods are the same
pR2(m3)
pR2(for.aic)
pR2(for.bic) # no diff really in r2

anova(m3, for.aic, for.bic)
regTermTest(m3, "ekm") # could well lose him

m4<-glm(PA~shd+chl_log+wnd+poly(tmc,2)+
          smt_sqt+poly(bty,2)+bet_adu_03_ave+bet_juv_03_ave+
          poly(yft_adu_03_ave,2)+skj_juv_03_ave+RAC,
        data=dat_heron, family="binomial")


# looking at variable importance/contribution
library(caret)
varImp(m4) 
summary(m4)
# all looks good but is showing a negative trend rather than
# a positive so is modelled wrong, this is probably due to collinearity
# make sure it is not modelling correct


# Check collinearity..need to see whats acceptable..could potentially remove 
# shd, tmc etc but pretty ruthless


temp<-dat_heron[1,-4]
temp[1,]<-as.vector(apply(dat_heron[,-4], 2,median))

heron_pred<-data.frame(yft_adu_03_ave=dat_heron$yft_adu_03_ave, 
                       temp)

p1<-predict(m3, newdata=heron_pred, type="response")

d1<-data.frame(PA=dat_heron$PA, yft_adu_03_ave=heron_pred$yft_adu_03_ave, pred=p1)
g1<-ggplot(data=d1, aes(y=PA, x=yft_adu_03_ave))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=pred, x=yft_adu_03_ave), colour=2)

# yeh looks bad, see if removing it makes much difference
m4<-update(m3, ~.- bet_juv_03_ave)
anova(m3, m4)
pR2(m3);pR2(m4) # doesnt add much at all

summary(m4) # now the tmc poly looks insignificant

m5<-glm(PA~chl_log+wnd+tmc+smt_sqt+poly(bty, 2)+
          bet_adu_03_ave+
          poly(skj_juv_03_ave,2),
        data=dat_heron, family="binomial")

anova(m4, m5)
pR2(m4);pR2(m5) # yep get rid, go with m5

# ok I think we're there
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
print(sum((resid(m5, type="pearson")^2))/df.residual(m5))
library(verification)
print(roc.area(dat_heron$PA, fitted(m1))$A)
print(roc.area(dat_heron$PA, fitted(m5))$A)

summary(m5)
varImp(m5) # all looks good, could try removing tmc all together?

m6<-update(m5, ~.- tmc)
anova(m5, m6)
pR2(m5);pR2(m6) # doesnt add much but
regTermTest(m5, "tmc") # is still significant so keep

library(hier.part) # see if we can kinda validate with hier.part
hp<-hier.part(y=dat_heron$PA, x=dat_heron[,c(8:11,15, 20, 21) ], family="binomial")
# ok so kinda the same but hier.part is sometimes dodgy:
#http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011698

m5_interp<-glm(PA~chl_log+wnd+tmc+smt_sqt+poly(bty, 2, raw=T)+
                 bet_adu_03_ave+
                 poly(skj_juv_03_ave,2, raw=T),
               data=dat_heron, family="binomial")
anova(m5, m5_interp)

summary(m5_interp)

