# 25/09/16 Lancaster, UK
# Model climatology data comparing foraging areas against background 

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

dat<-read.csv("spreads/paper2_extractionV4.csv", h=T, strip.white=T)
# V3 takes all PA points

d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "dtyp"))

p<-ggplot(d1, aes(dset, value))
p+geom_boxplot()+facet_wrap(~variable, scale="free")

## We're going to use the Information Theoretic approach to test presence (50%UD) absence (max range)
## of WTSH using a glm testing oceanographic vs tuna variables using AIC. We will do this 
## seperately for each colony.

# first we need to remove outliers and transform some variables

qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")

qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

qplot(Longitude, Latitude, data=dat)+geom_point(data=dat[dat$chl<1500,], aes(Longitude, Latitude), colour=2)

dat[!is.na(dat$chl) & dat$chl>1500,]$chl<-NA
# also knock out some dodgy looking wind values
dat[is.na(dat$chl),]$wnd<-NA
# and some weird sst
dat[!is.na(dat$sst) & dat$sst<0,]$sst<-NA
#weird ekem point
dat[!is.na(dat$ekm) &dat$ekm< -0.00001,]$ekm<-NA
#remove skj_juv outliers
dat[which(dat$skj_juv_03_ave>28),]$skj_juv_03_ave<-NA

## NOW we do naomit! 
## This is integral to modelling with the IT approach, because
## when subset models are created, they have different numbers of NAs
## and therefore data: they are not comparable

nrow(dat)
nrow(na.omit(dat))
dat<-na.omit(dat)
# mop up some presence stragglers that get cut by the na.omit
dat<-dat[-which(dat$dset=="LTHeronPTT2013" & dat$Longitude<152.61 &
                  dat$Latitude<(-18.76)),]

dat<-dat[-which(dat$dset=="LTHeronPTT2013" & dat$Longitude<151.62),]

dat<-dat[-which(dat$dset=="LTHeronGPS2015" &
                  dat$Longitude<153.24),]

# how we looking now
d1<-melt(dat, id.vars=c("Latitude", "Longitude", "dset", "dtyp"))
qplot(x=value, data=d1)+facet_wrap(~variable, scale="free")
qplot(factor(variable), value, data=d1, geom="boxplot")+facet_wrap(~variable, scale="free")

#g<-ggplot(data=d1, aes(x=factor(variable), y=value))
#g+geom_violin(alpha=0.5, color="gray")+
#  geom_jitter(alpha=0.5,position = position_jitter(width = 0.1))+
#  facet_wrap(~variable, scale="free")


# seamounts and chl still need some attention
qplot(chl, data=dat);qplot(log(chl), data=dat)
qplot(smt, data=dat);qplot(sqrt(smt), data=dat)

dat$chl_log<-log(dat$chl)
dat$smt_sqt<-sqrt(dat$smt)
dat$PA<-0
dat[dat$dtyp=="ud50_pres",]$PA<-1 # making presence absence variable

# setup colony datasets
# here we randomly assign years of tracking data to psuedo absences
# at a 1:10 ratio

####### Heron Island analyses #######

heronPA<-dat[dat$dset=="heronBuff",]

n=10 

dat_heron<-rbind(dat[dat$dset=="LTHeronGPS2015",], 
                 heronPA[sample(1:nrow(heronPA), 
                                (nrow(dat[dat$dset=="LTHeronGPS2015",])*n)),],
                 dat[dat$dset=="LTHeronPTT2011",], 
                 heronPA[sample(1:nrow(heronPA), 
                                (nrow(dat[dat$dset=="LTHeronPTT2011",])*n)),],
                 dat[dat$dset=="LTHeronPTT2013",], 
                 heronPA[sample(1:nrow(heronPA), 
                                (nrow(dat[dat$dset=="LTHeronPTT2013",])*n)),])

dat_heron$YRID<-factor(c(rep(2015, (nrow(dat_heron[dat_heron$dset==
                                                     "LTHeronGPS2015",])*(n+1))),
                         rep(2011, (nrow(dat_heron[dat_heron$dset==
                                                     "LTHeronPTT2011",])*(n+1))),
                         rep(2013, (nrow(dat_heron[dat_heron$dset==
                                                     "LTHeronPTT2013",])*(n+1)))))

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


m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          bet_juv_03_ave+yft_adu_03_ave+skj_juv_03_ave+YRID,
        data=dat_heron, family="binomial")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_heron$PA, fitted(m1))$A)

# now we see which variables are better suited to a second deg polynomial

d1<-melt(dat_heron, id.vars=c("PA", "YRID"))
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

# global model

m1<-glm(PA~poly(shd,2)+ekm+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          bet_adu_03_ave+bet_juv_03_ave+poly(yft_adu_03_ave,2)+
          poly(skj_juv_03_ave,2)+YRID,
        data=dat_heron, family="binomial")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_heron$PA, fitted(m1))$A)
library(pscl)
pR2(m1)

## spatial autocorrelation testing ##
library(pgirmess)
pgir1<-correlog(coords=dat_heron[,1:2], z=fitted(m1), method="Moran")

library(sp)
dh_ndub<-NULL
for(j in unique(dat_heron$YRID))
{ d1<-dat_heron[dat_heron$YRID==j,]
pts <- SpatialPoints(d1[,1:2])
pts <- SpatialPointsDataFrame(pts, data=d1)
dh_ndub<-rbind(dh_ndub, remove.duplicates(pts)@data)
print(j)}

nrow(dat_heron); nrow(dh_ndub)
## All points

# ok we know there is some multicollinearity in this model but it 
# gets worse as we select the model
library(gee)

gee1<-gee(PA~poly(shd,2)+ekm+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
            bet_adu_03_ave+bet_juv_03_ave+poly(yft_adu_03_ave,2)+
            poly(skj_juv_03_ave,2), id=YRID,
          data=dat_heron, family="binomial", corstr="independence")
## the above dies!
library(MASS)
library(nlme)

dh_ndub2<-dh_ndub[sample(1:nrow(dh_ndub), 1000),]

attach(dh_ndub2)
pql1<-glmmPQL(PA~poly(shd,2)+ekm+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
            bet_adu_03_ave+bet_juv_03_ave+poly(yft_adu_03_ave,2)+
            poly(skj_juv_03_ave,2), random=~1|YRID,
            data=dh_ndub2, family="binomial",
            correlation=corExp(form=~Longitude+Latitude))

detach(dh_ndub2)

# looking at variable importance/contribution
summary(m1) # look to drop ekm and yrid and probs yft
anova(m1) # anova tests terms sequentially
library(survey)
regTermTest(m1, "ekm")
anova(m1, update(m1, ~.- ekm)) # basically tells the same as anova but gives significance too

library(car)
Anova(m1) # Anova from car tests terms according to marginality
# i.e. after all other terms have been included
# not much support from ekm, YRID or yft

#### DAnGER !!!

dat_heron<-dat_heron[sample(1:nrow(dat_heron), 1000),]

m2<-glm(PA~poly(shd,2)+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          bet_adu_03_ave+bet_juv_03_ave+
          poly(skj_juv_03_ave,2),
        data=dat_heron, family="binomial")

# ok lets try gam and autocoavariate term models to incorp SPAC
library(ncf)
library(spdep)

RAC<-autocov_dist(resid(m2, type="pearson"), cbind(dat_heron[,1], dat_heron[,2]),
              nbs = 100, type = "inverse", zero.policy = T,
             style = "B", longlat=TRUE)

dat_heron<-cbind(dat_heron, RAC)

m3<-glm(PA~poly(shd,2)+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          bet_adu_03_ave+bet_juv_03_ave+
          poly(skj_juv_03_ave,2)+RAC,
        data=dat_heron, family="binomial") # throws 0 or 1 error


# so try gam
m3<-gam(PA~poly(shd,2)+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          bet_adu_03_ave+bet_juv_03_ave+
          poly(skj_juv_03_ave,2)+s(Longitude, Latitude),
        data=dat_heron, family="binomial") # throws 0 or 1 error

glm2corr<-spline.correlog(x=dat_heron$Longitude, y=dat_heron$Latitude, z=resid(m2, type="pearson"), latlon=T, resamp=1)

gam3corr<-spline.correlog(x=dat_heron$Longitude, y=dat_heron$Latitude, z=resid(m3, type="pearson"), latlon=T, resamp=1)



anova(m1, m2); AIC(m1, m2)

Anova(m2) # all appear significant

# See if step methods give same result from the m1 model

for.aic <- step(glm(PA~1,data=dat_heron, family="binomial"),
                direction = "forward", scope = formula(m1), k = 2, trace = 1) # forward AIC
for.bic <- step(glm(PA~1,data=dat_heron, family="binomial"),
                direction = "forward", scope = formula(m1), k = log(nrow(dat_heron)), trace = 0) # forward BIC
back.aic <- step(m1, direction = "backward", k = 2, trace = 0) # backward AIC
back.bic <- step(m1, direction = "backward", k = log(nrow(dat_heron)), trace = 0) # backward BIC

formula(for.aic);formula(for.bic);formula(back.aic);formula(back.bic);
formula(m2)
# so basically m2 and aic methods are the same, bic says remove shd also
pR2(m2)
pR2(for.aic)
pR2(for.bic) # no diff really in r2

anova(m2, for.aic, for.bic)
library(survey)
regTermTest(m2, "poly(shd, 2)") # could well lose him

m3<-glm(PA~chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          bet_adu_03_ave+bet_juv_03_ave+
          poly(skj_juv_03_ave,2),
        data=dat_heron, family="binomial")


# looking at variable importance/contribution
library(caret)
varImp(m3) 
summary(m3)
# all looks good but bet_ju_03_ave is showing a negative trend rather than
# a positive so is modelled wrong, this is probably due to collinearity
# make sure it is not modelling correct

temp<-dat_heron[1,-4]
temp[1,]<-as.vector(apply(dat_heron[,-4], 2,median))

heron_pred<-data.frame(bet_juv_03_ave=dat_heron$bet_juv_03_ave, 
                       temp,YRID="2013")

p1<-predict(m3, newdata=heron_pred, type="response")

d1<-data.frame(PA=dat_heron$PA, bet_juv_03_ave=heron_pred$bet_juv_03_ave, pred=p1)
g1<-ggplot(data=d1, aes(y=PA, x=bet_juv_03_ave))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=pred, x=bet_juv_03_ave), colour=2)

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

## attempt to use Dominance Analysis to ascertain variable importance
library(devtools)
install_github("clbustos/dominanceAnalysis")
library(dominanceanalysis)
# load in the below functions manually as dont seem to load in package 

dominanceMatrix<-function(x,undefined.value=0.5) {
  vars<-colnames(x)
  m<-length(vars)
  ma<-matrix(undefined.value,m,m,dimnames=list(vars,vars))
  
  for(i in 1:(m-1)) {
    for(j in (i+1):m) {
      comps<-na.omit(cbind(x[,i,drop=F],x[,j,drop=F]))
      if(mean(comps[,1]>comps[,2])==1) 
      {
        ma[i,j]<-1
        ma[j,i]<-0
      }
      
      if(mean(comps[,1]<comps[,2])==1) 
      {
        ma[i,j]<-0
        ma[j,i]<-1
      }
      
    }
  }
  ma
}

daAverageContributionByLevel<-function(x) {
  ff<-x$fit.functions
  
  out<-list()
  for(i in ff) {
    
    res<-aggregate(x$fits[[i]],list(level=x$level),mean,na.rm=T)
    
    out[[i]]<-res[res$level<max(x$level),]
    
  }
  out
}

daCompleteDominance<-function(daRR) {
  lapply(daRR$fits,dominanceMatrix)
}

daConditionalDominance<-function(daRR) {
  daACBL<-daAverageContributionByLevel(daRR)
  analize<-function(x) {
    x<-x[,-1]
    dominanceMatrix(x)
  }
  lapply(daACBL,analize)
}

daGeneralDominance<-function(daRR) {
  daACBL<-daAverageContributionByLevel(daRR)
  analize<-function(x) {
    x<-x[,-1]
    gm<-matrix(colMeans(x),1,ncol(x),dimnames=list(1,colnames(x)))
    dominanceMatrix(gm)
  }
  lapply(daACBL,analize)
}

dm<-dominanceAnalysis(m5)# fails

#slightly hacked version, had to remove 2 of the r2
# metrics (Nagelkerke and Cox and Snell) as throwing NAs

dm2<-function (x, constants = c(), fit.functions = "default", data = NULL, 
          null.model = NULL, ...) 
{
  daModels <- daSubmodels(x, constants)
  daRaw4r2 <- daRawResults(x, constants, fit.functions, data, 
                        null.model, ...)
  daRaw<-daRaw4r2
  daRaw$fit.functions<-daRaw$fit.functions[c(1,4)]
  daRaw$fits<-daRaw$fits[c(1,4)]
  daRaw$base.fits<-daRaw$base.fits[,c(1,4)]
  
  daAverageByLevel <- daAverageContributionByLevel(daRaw)
  daAverageGeneral <- lapply(daAverageByLevel, function(x) {
    colMeans(x[, -1])
  })
  list(predictors = daModels$predictors, constants = daModels$constants, 
       fit.functions = daRaw$fit.functions, fits = daRaw, contribution.by.level = daAverageByLevel, 
       contribution.average = daAverageGeneral, complete = daCompleteDominance(daRaw), 
       conditional = daConditionalDominance(daRaw), general = daGeneralDominance(daRaw))
}

dm<-dm2(m5) # works!
#I'll be using McFaddens' r2
dm$contribution.average$r2.m

#make sure it lines up
sum(dm$contribution.average$r2.m)
pR2(m5) # yep :)

#library(relaimpo)
#calc.relimp(m2)  only work with guassian link

# ok as I understand it the importance of the variables can be determined from the
#coefficients IF they are all scaled, and the model is not doing anything funky.
# to further investigate I need to keep googling logistic regression!

# Information Theoretic approach

m_intercept<-glm(PA~1, data=dat_heron, family="binomial")

m_oceano<-glm(PA~poly(shd,2)+ekm+chl_log+wnd+poly(tmc, 2)+smt_sqt+
                poly(bty,2), data=dat_heron, family="binomial")

m_tuna<-glm(PA~poly(bet_adu_03_ave,2)+bet_juv_03_ave+
              poly(yft_adu_03_ave,2)+poly(skj_juv_03_ave,2),
            data=dat_heron, family="binomial")

m_shd<-glm(PA~poly(shd, 2), data=dat_heron, family="binomial")
m_ekm<-glm(PA~ekm,data=dat_heron, family="binomial")
m_chl<-glm(PA~chl_log, data=dat_heron, family="binomial")
m_wnd<-glm(PA~wnd, data=dat_heron, family="binomial")
m_tmc<-glm(PA~poly(tmc, 2), data=dat_heron, family="binomial")
m_smt<-glm(PA~smt_sqt, data=dat_heron, family="binomial")
m_bty<-glm(PA~poly(bty, 2), data=dat_heron, family="binomial")
m_betadu<-glm(PA~poly(bet_adu_03_ave, 2),data=dat_heron, family="binomial")
m_betjuv<-glm(PA~bet_juv_03_ave, data=dat_heron, family="binomial")
m_yft<-glm(PA~poly(yft_adu_03_ave, 2), data=dat_heron, family="binomial")
m_skjjuv<-glm(PA~poly(skj_juv_03_ave,2),data=dat_heron, family="binomial")
m_yr<-glm(PA~YRID,data=dat_heron, family="binomial")

temp<-AIC(m_intercept, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
          m_betadu, m_betjuv, m_yft, m_skjjuv, m_yr)

mod.list<-list(m_intercept, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
               m_betadu, m_betjuv, m_yft, m_skjjuv, m_yr)

it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out

# additional trial with biologically plausable models
m_intercept<-glm(PA~1, data=dat_heron, family="binomial")
m_topo<-glm(PA~smt_sqt+poly(bty, 2), data=dat_heron, family="binomial")
m_oceo<-glm(PA~poly(tmc, 2)+poly(shd, 2), data=dat_heron, family="binomial")
m_prod<-glm(PA~chl_log, data=dat_heron, family="binomial")
m_wnd<-glm(PA~wnd, data=dat_heron, family="binomial")
m_yr<-glm(PA~YRID,data=dat_heron, family="binomial")
m_liltun<-glm(PA~bet_juv_03_ave+poly(skj_juv_03_ave,2),
            data=dat_heron, family="binomial")
m_bigtun<-glm(PA~poly(bet_adu_03_ave,2)+  poly(yft_adu_03_ave,2),
              data=dat_heron, family="binomial")

temp<-AIC(m_intercept, m_topo, m_oceo, m_prod, m_wnd,  
           m_liltun, m_bigtun, m_yr)

mod.list<-list(m_intercept, m_topo, m_oceo, m_prod, m_wnd,  
               m_liltun, m_bigtun, m_yr)

it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out

# hmm try the dredge approach

options(na.action = na.fail)
GD1 <- dredge(m1, m.lim = c(NA,3))
# I then write out subset (GD1, delta<300) for reference

#######$ NOW for Lord Howe analyses ########

lhiPA<-dat[dat$dset=="lhiBuff",]
 
n=10

dat_lhi<-rbind(dat[dat$dset=="LTLHIGPS2014",], 
                 lhiPA[sample(1:nrow(lhiPA), 
                                (nrow(dat[dat$dset=="LTLHIGPS2014",])*n)),],
                 dat[dat$dset=="LTLHIGPS2015",], 
                 lhiPA[sample(1:nrow(lhiPA), 
                                (nrow(dat[dat$dset=="LTLHIGPS2015",])*n)),],
                 dat[dat$dset=="LTLHIGPS2016",], 
                 lhiPA[sample(1:nrow(lhiPA), 
                                (nrow(dat[dat$dset=="LTLHIGPS2016",])*n)),])

dat_lhi$YRID<-factor(c(rep(2014, (nrow(dat_lhi[dat_lhi$dset==
                                                     "LTLHIGPS2014",])*(n+1))),
                         rep(2015, (nrow(dat_lhi[dat_lhi$dset==
                                                     "LTLHIGPS2015",])*(n+1))),
                         rep(2016, (nrow(dat_lhi[dat_lhi$dset==
                                                     "LTLHIGPS2016",])*(n+1)))))

# z-transform (scale and center) to make variables comparable on same scale
# nope
#library(vegan)
#dat_lhi<-cbind(dat_lhi[,c(1,2,24:25)],
#                 decostand(dat_lhi[,c(5:7, 9,10,12:23)], method="standardize"))

# now have a look at collinearity
library(car)
vif(lm(1:nrow(dat_lhi)~sst+shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+bet_tot_03_ave+yft_adu_03_ave+
         yft_juv_03_ave+yft_tot_03_ave+skj_adu_03_ave+
         skj_juv_03_ave+skj_tot_03_ave, data=dat_lhi))  

vif(lm(1:nrow(dat_lhi)~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
         bet_juv_03_ave+yft_adu_03_ave+
         skj_adu_03_ave,
       data=dat_lhi))

# ok looks good, all values < 10, hmm skjadu and yftadu are corr according to 
cor(dat_lhi[,-(1:4)])

#pairs(dat_lhi[,-(1:2)], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_lhi[,11:19], upper.panel = panel.smooth,lower.panel=panel.cor)
# ok so actually for lhi I include adult yft and bet and juv skj, 
# juv skj are proxy for adult skj, juv yft and juv bet.

m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          yft_adu_03_ave+skj_juv_03_ave+YRID,
        data=dat_lhi, family="binomial")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_lhi$PA, fitted(m1))$A)

# now we see which variables are better suited to a second deg polynomial

d1<-melt(dat_lhi, id.vars=c("PA", "YRID"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=PA, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")

# we would do below for each variable individually to see the suitability of a poly

m_lin<-glm(PA~skj_juv_03_ave,data=dat_lhi, family="binomial")
m_pol<-glm(PA~poly(skj_juv_03_ave, 2),data=dat_lhi, family="binomial")
d1<-data.frame(PA=dat_lhi$PA, skj_juv_03_ave=dat_lhi$skj_juv_03_ave, m_lin=fitted(m_lin), m_pol=fitted(m_pol))
g1<-ggplot(data=d1, aes(y=PA, x=skj_juv_03_ave))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=m_lin, x=skj_juv_03_ave), colour=2)+
  geom_line(aes(y=m_pol, x=skj_juv_03_ave), colour=3)
anova(m_lin, m_pol)
#poly model better

# global model

m1<-glm(PA~poly(shd,2)+ekm+poly(chl_log,2)+poly(wnd,2)+tmc+
          poly(smt_sqt,2)+poly(bty, 2)+bet_adu_03_ave+
          poly(yft_adu_03_ave,2)+poly(skj_juv_03_ave,2)+YRID,
        data=dat_lhi, family="binomial")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_lhi$PA, fitted(m1))$A)

## spatial autocorrelation testing ##
library(pgirmess)
pgir1<-correlog(coords=dat_lhi[,1:2], z=fitted(m1), method="Moran")

# we notied that in m1 polys of shd, tmc, wnd and skj and bet threw
# the 0 or 1 error. bet and tmc didnt need a poly, skj is problematic but I think
# it needs the poly. shd and wnd could use polys but wont work so meh
m2<-glm(PA~shd+ekm+poly(chl_log,2)+wnd+tmc+
          poly(smt_sqt,2)+poly(bty, 2)+bet_adu_03_ave+
          yft_adu_03_ave+poly(skj_juv_03_ave,2)+YRID,
        data=dat_lhi, family="binomial")

anova(m1,m2) # ok non poly model worse but not by anything much

pairs(dat_lhi[,c(6:11,15,17, 20, 21) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# looking at it in full, i could be tempted to remove shd and yellowfin or BET adult 

# looking at variable importance/contribution

library(car)
Anova(m2) # Anova from car tests terms according to marginality
# i.e. after all other terms have been included
# not much support from wnd or YRID.. but beware the error
m3<-glm(PA~shd+ekm+poly(chl_log,2)+tmc+
          poly(smt_sqt,2)+poly(bty, 2)+bet_adu_03_ave+
          yft_adu_03_ave+poly(skj_juv_03_ave,2),
        data=dat_lhi, family="binomial")

anova(m1, m2,m3); AIC(m1, m2, m3)

Anova(m2) # all appear significant

# See if step methods give same result from the m1 model

for.aic <- step(glm(PA~1,data=dat_lhi, family="binomial"),
                direction = "forward", scope = formula(m2), k = 2, trace = 1) # forward AIC
for.bic <- step(glm(PA~1,data=dat_lhi, family="binomial"),
                direction = "forward", scope = formula(m2), k = log(nrow(dat_lhi)), trace = 0) # forward BIC
back.aic <- step(m2, direction = "backward", k = 2, trace = 0) # backward AIC
back.bic <- step(m2, direction = "backward", k = log(nrow(dat_lhi)), trace = 0) # backward BIC

formula(for.aic);formula(for.bic);formula(back.aic);formula(back.bic);
formula(m3)
# so basically m3 and bic methods are the same, aic says leave YRID and wnd
pR2(m2)
pR2(m3)
pR2(for.aic)
pR2(for.bic) # no diff really in r2

anova(m2, m3, for.aic, for.bic)
library(survey)
regTermTest(m2, "wnd") # could well lose him
regTermTest(m2, "YRID") # could well lose him

# m3 it is!

# looking at variable importance/contribution
library(caret)
varImp(m3) 
summary(m3)
# all looks good could probably lose that chl poly too

m4<-glm(PA~shd+ekm+chl_log+tmc+
          poly(smt_sqt,2)+poly(bty, 2)+bet_adu_03_ave+
          yft_adu_03_ave+poly(skj_juv_03_ave,2),
        data=dat_lhi, family="binomial")

anova(m3, m4); AIC(m3, m4);pR2(m3);pR2(m4) # yeh I'll lose him

# ok I think we're there
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
print(sum((resid(m4, type="pearson")^2))/df.residual(m4))
library(verification)
print(roc.area(dat_lhi$PA, fitted(m1))$A)
print(roc.area(dat_lhi$PA, fitted(m4))$A)

summary(m4)

m4_interp<-glm(PA~shd+ekm+chl_log+tmc+
                 poly(smt_sqt,2, raw=T)+poly(bty, 2, raw=T)+bet_adu_03_ave+
                 yft_adu_03_ave+poly(skj_juv_03_ave,2, raw=T),
               data=dat_lhi, family="binomial")

summary(m4_interp)

# using my hacked 'dm2' dominanceAnalysis function from heron analyses above
library(dominanceanalysis)
dm<-dm2(m4) # works!
#I'll be using McFaddens' r2
dm$contribution.average$r2.m

#make sure it lines up
sum(dm$contribution.average$r2.m)
pR2(m4) # yep :)

# I think the skj poly in dm2 gives the familiar glm.fit 1-0 warnings, 
# just to confirm to myself here is same model withou poly
m5<-glm(PA~shd+ekm+chl_log+tmc+
          poly(smt_sqt,2)+poly(bty, 2)+bet_adu_03_ave+
          yft_adu_03_ave+skj_juv_03_ave,
        data=dat_lhi, family="binomial")

pR2(m4);pR2(m5) #yeh ok its worse we know
dm5<-dm2(m5) # still get 1 message! but not like 50 ;)

dm5$contribution.average$r2.m
dm$contribution.average$r2.m



# Information Theoretic approach

m_intercept<-glm(PA~1, data=dat_lhi, family="binomial")

m_oceano<-glm(PA~poly(shd,2)+ekm+poly(chl_log,2)+poly(wnd,2)+
              poly(tmc,2)+poly(smt_sqt,2)+poly(bty,2),
              data=dat_lhi, family="binomial")

m_tuna<-glm(PA~poly(bet_adu_03_ave,2)+
              yft_adu_03_ave+poly(skj_juv_03_ave,2),
            data=dat_lhi, family="binomial")

m_shd<-glm(PA~poly(shd, 2), data=dat_lhi, family="binomial")
m_ekm<-glm(PA~ekm,data=dat_lhi, family="binomial")
m_chl<-glm(PA~poly(chl_log,2), data=dat_lhi, family="binomial")
m_wnd<-glm(PA~poly(wnd,2), data=dat_lhi, family="binomial")
m_tmc<-glm(PA~poly(tmc, 2), data=dat_lhi, family="binomial")
m_smt<-glm(PA~poly(smt_sqt,2), data=dat_lhi, family="binomial")
m_bty<-glm(PA~poly(bty, 2), data=dat_lhi, family="binomial")
m_betadu<-glm(PA~poly(bet_adu_03_ave, 2),data=dat_lhi, family="binomial")
m_yft<-glm(PA~yft_adu_03_ave, data=dat_lhi, family="binomial")
m_skjjuv<-glm(PA~poly(skj_juv_03_ave,2),data=dat_lhi, family="binomial")
m_yr<-glm(PA~YRID,data=dat_lhi, family="binomial")

temp<-AIC(m_intercept, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
          m_betadu, m_yft, m_skjjuv, m_yr)

mod.list<-list(m_intercept, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
               m_betadu, m_yft, m_skjjuv, m_yr)

it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out

## playing with plotting model results

m3<-glm(PA~shd+ekm+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          poly(bet_adu_03_ave,2)+bet_juv_03_ave+poly(yft_adu_03_ave,2)+
          poly(skj_juv_03_ave,2)+YRID,
        data=dat_heron, family="binomial")

temp<-dat_heron[1,-4]
temp[1,]<-as.vector(apply(dat_heron[,-4], 2,median))

heron_pred<-data.frame(shd=dat_heron$shd, 
                       temp,YRID="2013")

p1<-predict(m1, newdata=heron_pred, type="response")

d1<-data.frame(PA=dat_heron$PA, shd=heron_pred$shd, pred=p1)
g1<-ggplot(data=d1, aes(y=PA, x=shd))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=pred, x=shd), colour=2)

#poly model better


##### BELOW IS OLD GAM STUFF!!! ### kept in here for reference

library(mgcv)
m2<-gam(PA~s(shd, k=3)+s(ekm, k=3)+s(chl_log, k=3)+s(wnd, k=3)+s(tmc, k=3)+
          s(smt_sqt, k=3)+s(bty, k=3)+s(bet_adu_03_ave, k=3)+
          s(bet_juv_03_ave, k=3)+s(yft_adu_03_ave, k=3)+s(skj_juv_03_ave, k=3)+YRID,
        data=dat_heron, family="binomial")
summary(m2)
print(sum((resid(m2, type="pearson")^2))/df.residual(m2))
library(verification)
print(roc.area(dat_heron$PA, fitted(m2))$A)

m_intercept<-gam(PA~1, data=na.omit(dat_heron), family="binomial")

m_oceano<-gam(PA~s(shd, k=3)+s(ekm, k=3)+s(chl_log, k=3)+s(wnd, k=3)+s(tmc, k=3)+s(smt_sqt, k=3)+s(bty, k=3),
              data=dat_heron, family="binomial")

m_tuna<-gam(PA~s(bet_adu_03_ave, k=3)+s(bet_juv_03_ave, k=3)+s(yft_adu_03_ave,k=3)+s(skj_juv_03_ave, k=3),
            data=dat_heron, family="binomial")

m_shd<-gam(PA~s(shd, k=3), data=dat_heron, family="binomial")
m_ekm<-gam(PA~s(ekm, k=3),data=dat_heron, family="binomial")
m_chl<-gam(PA~s(chl_log, k=3), data=dat_heron, family="binomial")
m_wnd<-gam(PA~s(wnd, k=3), data=dat_heron, family="binomial")
m_tmc<-gam(PA~s(tmc, k=3), data=dat_heron, family="binomial")
m_smt<-gam(PA~s(smt_sqt, k=3), data=dat_heron, family="binomial")
m_bty<-gam(PA~s(bty, k=3), data=dat_heron, family="binomial")
m_betadu<-gam(PA~s(bet_adu_03_ave, k=3),data=dat_heron, family="binomial")
m_betjuv<-gam(PA~s(bet_juv_03_ave, k=3), data=dat_heron, family="binomial")
m_yft<-gam(PA~s(yft_adu_03_ave, k=3), data=dat_heron, family="binomial")
m_skjjuv<-gam(PA~s(skj_juv_03_ave, k=3),data=dat_heron, family="binomial")
m_mo<-gam(PA~YRID,data=dat_heron, family="binomial")


temp<-AIC(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
          m_betadu, m_yft, m_skjjuv, m_mo)


mod.list<-list(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
               m_betadu, m_yft, m_skjjuv, m_mo)


it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out







d1<-melt(dat_heron, id.vars=c("PA", "YRID"))

g1<-ggplot(data=d1, aes(y=PA, x=value))
g1+geom_jitter(height=0.1)+geom_smooth(method="glm")+facet_wrap(~variable, scale="free")

g1+geom_jitter(aes(colour=YRID),height=0.1, size=0.5)+geom_smooth(formula=y~s(x, k=5),method="gam")+facet_wrap(~variable, scale="free")
# 

dat[dat$dset=="heronBuff"  | dat$dset=="LTHeronGPS2015"  | 
      dat$dset=="LTHeronPTT2011"  | dat$dset=="LTHeronPTT2013" , ]

dat_lhi<-dat[dat$dset=="lhiBuff"  | dat$dset=="LTLHIGPS2014"  | 
               dat$dset=="LTLHIGPS2015"  | dat$dset=="LTLHIGPS2016" , ]

he1<-melt(dat_heron, id.vars=c("Latitude", "Longitude", "dset", "dtyp", "PA"))
qplot(factor(variable), value, colour= factor(PA), data=he1, geom="boxplot")+facet_wrap(~variable, scale="free")

lh1<-melt(dat_lhi, id.vars=c("Latitude", "Longitude", "dset", "dtyp", "PA"))
qplot(factor(variable), value, colour= factor(PA), data=lh1, geom="boxplot")+facet_wrap(~variable, scale="free")


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




#library(PerformanceAnalytics)
#chart.Correlation(dat_heron[,c(6,7,9,10,12,13,14,16,18,22)], histogram=F, pch=19)

mcor<-cor(na.omit(dat_heron[,-c(1:4)]))
mcor # we can see SKJ adu (0.8)and YFT adu (0.67) are highly coorelated with their juvs
mcor<-cor(na.omit(dat_lhi[,-c(1:4)]))
mcor #we can see SKJ adu (0.8) is highly coorelated with its juv, YFT less so (0.2)

### try some models

m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          bet_juv_03_ave+yft_adu_03_ave+skj_juv_03_ave+YRID,
        data=dat_heron, family="binomial")
summary(m1)
sum((resid(m1, type="pearson")^2))/df.residual(m1)

library(MuMIn)
options(na.action = "na.fail") 
dr<-dredge(m1)

#lets not jump the gun. just try some models first

m_intercept<-glm(PA~1, data=dat_heron, family="binomial")

m_tuna<-glm(PA~bet_adu_03_ave+bet_juv_03_ave+yft_adu_03_ave+skj_juv_03_ave,
            data=dat_heron, family="binomial")

m_oceano<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty,
              data=dat_heron, family="binomial")

m_shd<-glm(PA~shd, data=dat_heron, family="binomial")
m_ekm<-glm(PA~ekm,data=dat_heron, family="binomial")
m_chl<-glm(PA~chl_log, data=dat_heron, family="binomial")
m_wnd<-glm(PA~wnd, data=dat_heron, family="binomial")
m_tmc<-glm(PA~tmc, data=dat_heron, family="binomial")
m_smt<-glm(PA~smt_sqt, data=dat_heron, family="binomial")
m_bty<-glm(PA~bty, data=dat_heron, family="binomial")
m_betadu<-glm(PA~bet_adu_03_ave,data=dat_heron, family="binomial")
m_betjuv<-glm(PA~bet_juv_03_ave, data=dat_heron, family="binomial")
m_yft<-glm(PA~yft_adu_03_ave, data=dat_heron, family="binomial")
m_skj<-glm(PA~skj_juv_03_ave,data=dat_heron, family="binomial")

temp<-AIC(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
          m_betadu, m_betjuv, m_yft, m_skj)


mod.list<-list(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
               m_betadu, m_betjuv, m_yft, m_skj)


it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out



# If i run each year seperately the model becomes zero inflated. I can either sub
# sample the absence data to reduce the zinf then run individually or run with
# a random effect for dset. I guess i need to decide if i want to ask if the pattern
# varies between year, ie is it consistent or just the general pattern for each col... maybe the latter
# rather than ranef could include as factor.. the problem with both is allocating to the PA data, 
# needs to be stable.. can test

library(mgcv)

m1<-gam(PA~s(shd, k=3)+s(ekm, k=3)+s(chl_log, k=3)+s(wnd, k=3)+s(tmc, k=3)+s(smt_sqt, k=3)+s(bty, k=3)+
          s(bet_adu_03_ave, k=3)+s(bet_juv_03_ave, k=3)+s(yft_adu_03_ave,k=3)+s(skj_adu_03_ave, k=3)+YRID,
        data=dat_heron, family="binomial")
summary(m1)
sum((resid(m1, type="pearson")^2))/df.residual(m1)


m_intercept<-gam(PA~1, data=na.omit(dat_heron), family="binomial")

m_oceano<-gam(PA~s(shd, k=3)+s(ekm, k=3)+s(chl_log, k=3)+s(wnd, k=3)+s(tmc, k=3)+s(smt_sqt, k=3)+s(bty, k=3),
              data=na.omit(dat_heron), family="binomial")

m_tuna<-gam(PA~s(bet_adu_03_ave, k=3)+s(bet_juv_03_ave, k=3)+s(yft_adu_03_ave,k=3)+s(skj_adu_03_ave, k=3),
            data=na.omit(dat_heron), family="binomial")

m_shd<-gam(PA~s(shd, k=3), data=na.omit(dat_heron), family="binomial")
m_ekm<-gam(PA~s(ekm, k=3),data=na.omit(dat_heron), family="binomial")
m_chl<-gam(PA~s(chl_log, k=3), data=na.omit(dat_heron), family="binomial")
m_wnd<-gam(PA~s(wnd, k=3), data=na.omit(dat_heron), family="binomial")
m_tmc<-gam(PA~s(tmc, k=3), data=na.omit(dat_heron), family="binomial")
m_smt<-gam(PA~s(smt_sqt, k=3), data=na.omit(dat_heron), family="binomial")
m_bty<-gam(PA~s(bty, k=3), data=na.omit(dat_heron), family="binomial")
m_betadu<-gam(PA~s(bet_adu_03_ave, k=3),data=na.omit(dat_heron), family="binomial")
m_betjuv<-gam(PA~s(bet_juv_03_ave, k=3), data=na.omit(dat_heron), family="binomial")
m_yft<-gam(PA~s(yft_adu_03_ave, k=3), data=na.omit(dat_heron), family="binomial")
m_skj<-gam(PA~s(skj_adu_03_ave, k=3),data=na.omit(dat_heron), family="binomial")
m_skjjuv<-gam(PA~s(skj_juv_03_ave, k=3),data=na.omit(dat_heron), family="binomial")
m_yftjuv<-gam(PA~s(yft_juv_03_ave, k=3),data=na.omit(dat_heron), family="binomial")


temp<-AIC(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
          m_betadu, m_betjuv, m_yft, m_skj, m_skjjuv, m_yftjuv)


mod.list<-list(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
               m_betadu, m_betjuv, m_yft, m_skj, m_skjjuv, m_yftjuv)


it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out

m3<-gam(PA~s(shd, k=3)+s(ekm, k=3)+s(chl_log, k=3)+wnd+s(tmc, k=3)+s(smt_sqt, k=3)+s(bty, k=3)+
          bet_adu_03_ave+bet_juv_03_ave+yft_adu_03_ave+skj_juv_03_ave+YRID,
        data=dat_heron, family="binomial")

## ok so lets try modelling each year seperately while testing the effect of number of psuedo absences

### testing the number of psuedo absences

hb<-dat_heron[dat_heron$dset=="heronBuff",]
for(i in 1:50)
{
  d1<-rbind(dat_heron[dat_heron$dset=="LTHeronGPS2015",], hb[sample(1:nrow(hb), 20000),])
  m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
            bet_juv_03_ave+yft_adu_03_ave+skj_adu_03_ave,
          data=d1, family="binomial")
  print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
}
# ok so seems like we need more psedo absences and this will mean a zinf model..
install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")

summary(m1)
sum((resid(m1, type="pearson")^2))/df.residual(m1)


m1<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
          bet_juv_03_ave+yft_adu_03_ave+skj_adu_03_ave,
        data=dat_heron, family="binomial")

summary(m1)
sum((resid(m1, type="pearson")^2))/df.residual(m1)

m_intercept<-glm(PA~1, data=dat_heron, family="binomial")

m_tuna<-glm(PA~bet_adu_03_ave+bet_juv_03_ave+yft_adu_03_ave+skj_adu_03_ave,
            data=dat_heron, family="binomial")

m_oceano<-glm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty,
              data=dat_heron, family="binomial")

m_shd<-glm(PA~shd, data=dat_heron, family="binomial")
m_ekm<-glm(PA~ekm,data=dat_heron, family="binomial")
m_chl<-glm(PA~chl_log, data=dat_heron, family="binomial")
m_wnd<-glm(PA~wnd, data=dat_heron, family="binomial")
m_tmc<-glm(PA~tmc, data=dat_heron, family="binomial")
m_smt<-glm(PA~smt_sqt, data=dat_heron, family="binomial")
m_bty<-glm(PA~bty, data=dat_heron, family="binomial")
m_betadu<-glm(PA~bet_adu_03_ave,data=dat_heron, family="binomial")
m_betjuv<-glm(PA~bet_juv_03_ave, data=dat_heron, family="binomial")
m_yft<-glm(PA~yft_adu_03_ave, data=dat_heron, family="binomial")
m_skj<-glm(PA~skj_adu_03_ave,data=dat_heron, family="binomial")

temp<-AIC(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
          m_betadu, m_betjuv, m_yft, m_skj)


mod.list<-list(m_intercept, m_tuna, m_oceano, m_shd, m_ekm, m_chl, m_wnd, m_tmc, m_smt, m_bty, 
               m_betadu, m_betjuv, m_yft, m_skj)


it.out <- data.frame(loglik=sapply(mod.list, logLik), temp)
it.out$delta <- it.out$AIC-min(it.out$AIC)
it.out <- it.out[order(it.out$delta),]

it.out$w <- round(exp(-(1/2)*it.out$delta) / sum(exp(-(1/2)*it.out$delta)), 3)
it.out$cum.w <- cumsum(it.out$w)
it.out




#hmm ok so we're getting somewhere..
## try a brt for support. NOPE it kills R??
library(gbm)

heron.gbm<-gbm(PA~shd+ekm+chl_log+wnd+tmc+smt_sqt+bty+bet_adu_03_ave+
                 bet_juv_03_ave+yft_adu_03_ave+skj_adu_03_ave, data=na.omit(dat_heron),
               distribution="bernoulli", n.trees=10000,
               interaction.depth=3, train.fraction=0.75,
               bag.fraction=0.5, cv.folds=3, shrinkage=0.001,
               n.minobsinnode=2)
