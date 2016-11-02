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

dat<-read.csv("spreads/paper2_extractionV3.csv", h=T, strip.white=T)
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

library(vegan)
dat_heron<-cbind(dat_heron[,c(1,2,24:25)],
                 decostand(dat_heron[,c(5:7, 9,10,12:23)], method="standardize"))

# remove some other problematic datapoints
dat_heron<-dat_heron[-which(row.names(dat_heron)%in%c("72622", "72614", "72606")),]

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

m_lin<-glm(PA~shd,data=dat_heron, family="binomial")
m_pol<-glm(PA~poly(shd, 2),data=dat_heron, family="binomial")
d1<-data.frame(PA=dat_heron$PA, shd=dat_heron$shd, m_lin=fitted(m_lin), m_pol=fitted(m_pol))
g1<-ggplot(data=d1, aes(y=PA, x=shd))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=m_lin, x=shd), colour=2)+
  geom_line(aes(y=m_pol, x=shd), colour=3)
anova(m_lin, m_pol)
#poly model better

# global model

m1<-glm(PA~poly(shd,2)+ekm+chl_log+wnd+poly(tmc, 2)+smt_sqt+poly(bty, 2)+
          poly(bet_adu_03_ave,2)+bet_juv_03_ave+poly(yft_adu_03_ave,2)+
          poly(skj_juv_03_ave,2)+YRID,
        data=dat_heron, family="binomial")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_heron$PA, fitted(m1))$A)

## spatial autocorrelation testing ##
library(pgirmess)
pgir1<-correlog(coords=dat_heron[,1:2], z=fitted(m1), method="Moran")


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

library(vegan)
dat_lhi<-cbind(dat_lhi[,c(1,2,24:25)],
                 decostand(dat_lhi[,c(5:7, 9,10,12:23)], method="standardize"))

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

m_lin<-glm(PA~shd,data=dat_lhi, family="binomial")
m_pol<-glm(PA~poly(shd, 2),data=dat_lhi, family="binomial")
d1<-data.frame(PA=dat_lhi$PA, shd=dat_lhi$shd, m_lin=fitted(m_lin), m_pol=fitted(m_pol))
g1<-ggplot(data=d1, aes(y=PA, x=shd))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=m_lin, x=shd), colour=2)+
  geom_line(aes(y=m_pol, x=shd), colour=3)
anova(m_lin, m_pol)
#poly model better

# global model

m1<-glm(PA~poly(shd,2)+ekm+poly(chl_log,2)+poly(wnd,2)+poly(tmc, 2)+
          poly(smt_sqt,2)+poly(bty, 2)+poly(bet_adu_03_ave,2)+
          yft_adu_03_ave+poly(skj_juv_03_ave,2)+YRID,
        data=dat_lhi, family="binomial")
summary(m1)
print(sum((resid(m1, type="pearson")^2))/df.residual(m1))
library(verification)
print(roc.area(dat_lhi$PA, fitted(m1))$A)

## spatial autocorrelation testing ##
library(pgirmess)
pgir1<-correlog(coords=dat_lhi[,1:2], z=fitted(m1), method="Moran")

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
