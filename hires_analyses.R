# 10/10/16 Lancaster, UK
# Model HMM data modelling foraging against non-foraging


rm(list=ls())
library(ggplot2)
library(reshape2)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/GPS_LT_141516_hmm_oceano_attribs.csv", h=T, strip.white=T)

# remove a few pesky STs that sneaked thru
dat<-dat[dat$trip_id != "25_02_15_06LW3"& dat$trip_id != "28_03_15_04LW4"&
           dat$trip_id != "16_02_16_41LW3"& dat$trip_id !="18_02_16_22LW3",]

d1<-melt(dat[,c(1:5, 26:44)], id.vars=c("Latitude", "Longitude", "Date", "Time", "TrackID"))

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

# need to merge sst and chl datasets to try and reduce NAs
dat$SST<-apply(cbind(dat$sstAG,dat$sstOi), 1, FUN=mean,na.rm=T)
dat$ASST<-apply(cbind(dat$AsstAG,dat$AsstOi), 1, FUN=mean,na.rm=T)
dat$CHL_M<-apply(cbind(dat$chlMH1_M,dat$chlVH3_M), 1, FUN=mean,na.rm=T)
dat$CHL_A<-apply(cbind(dat$chlMH1,dat$chlVH3,dat$chlMH1_M,dat$chlVH3_M), 1, FUN=mean,na.rm=T)

nrow(dat)
for (i in 26:42)
{print(names(dat)[i]);print(length(which(is.na(dat[,i]))))}

#ok cool so not too many NAs if we use the right variables. 

#transformations
dat$LOGCHL_M<-log(dat$CHL_M)
dat$LOGCHL_A<-log(dat$CHL_A)
dat$SQRTsmt<-sqrt(dat$smt)

#collinearity
mcor<-cor(na.omit(dat[26:39]))
mcor 

library(car)
vif(lm(1:nrow(dat)~ekmU+modW+SST+ASST+
         sshHy+bty+SQRTsmt+LOGCHL_A+yftA+yftJ+skjA+skjJ+
         betA+betJ,data=dat)) # all good

#start modelling
# get HMM into binary response
dat$for_bin<-0
dat[dat$HMMstates!=3,]$for_bin<-1

#setup NA removed col year datasets
# Very important to remove NA prior to modelling
d_mod<-dat[,c(6,16,23,30,31,36:42,44:46,50:52)]
HER15<-na.omit(d_mod[d_mod$Colony=="Heron",])
LHI14<-na.omit(d_mod[d_mod$Colony=="LHI" & d_mod$Year==2014,])
LHI15<-na.omit(d_mod[d_mod$Colony=="LHI"& d_mod$Year==2015,])
LHI16<-na.omit(d_mod[d_mod$Colony=="LHI"& d_mod$Year==2016,])

vif(lm(1:nrow(HER15)~ekmU+modW+SST+ASST+
         sshHy+bty+SQRTsmt+LOGCHL_A+yftA+yftJ+skjA+skjJ+
         betA+betJ,data=HER15)) # all good


#removeunused factors
HER15$trip_id<-factor(HER15$trip_id)
LHI14$trip_id<-factor(LHI14$trip_id)
LHI15$trip_id<-factor(LHI15$trip_id)
LHI16$trip_id<-factor(LHI16$trip_id)

# do z transform
library(vegan)
HER15<-data.frame(HER15[,c(12,1,2,3)],decostand(HER15[,4:11],method="standardize"))
LHI14<-data.frame(LHI14[,c(12,1,2,3)],decostand(LHI14[,4:11],method="standardize"))
LHI15<-data.frame(LHI15[,c(12,1,2,3)],decostand(LHI15[,4:11],method="standardize"))
LHI16<-data.frame(LHI16[,c(12,1,2,3)],decostand(LHI16[,4:11],method="standardize"))


library(lme4) # initially with glmm
library(verificaton)
library(mgcv)
library(gamm4)
library(MuMIn)

d1<-melt(HER15, id.vars=c("for_bin","trip_id", "Colony", "Year"))

g1<-ggplot(data=d1, aes(y=for_bin, x=value))
g1+geom_jitter(height=0.1)+geom_smooth(method="glm")+facet_wrap(~variable, scale="free")
g1+geom_jitter(aes(colour=trip_id),height=0.1, size=0.5)+geom_smooth(formula=y~s(x, k=3),method="gam")+facet_wrap(~variable, scale="free")
# gam the right way to go, check k value
he15glm<-glm(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
              SQRTsmt+LOGCHL_A, family="binomial", 
            data=HER15,na.action=na.omit)


he15glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
            SQRTsmt+LOGCHL_A + (1|trip_id), family="binomial", 
            data=HER15,na.action=na.omit)

summary(he15glm)
sum((resid(he15glm, type="pearson")^2))/df.residual(he15glm)
roc.area(obs=HER15$for_bin,
         pred=fitted(he15glm)) 

summary(he15glmer)
sum((resid(he15glmer, type="pearson")^2))/df.residual(he15glmer)
roc.area(obs=HER15$for_bin,pred=fitted(he15glmer)) 

r.squaredGLMM(he15glmer)

lh14glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
                   SQRTsmt+LOGCHL_A + (1|trip_id), family="binomial", 
                 data=LHI14,na.action=na.omit)
summary(lh14glmer)
sum((resid(lh14glmer, type="pearson")^2))/df.residual(lh14glmer)

lh15glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
                   SQRTsmt+LOGCHL_A + (1|trip_id), family="binomial", 
                 data=LHI15,na.action=na.omit)
summary(lh15glmer)
sum((resid(lh15glmer, type="pearson")^2))/df.residual(lh15glmer)

lh16glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
                   SQRTsmt+LOGCHL_A + (1|trip_id), family="binomial", 
                 data=LHI16,na.action=na.omit)
summary(lh16glmer)
sum((resid(lh16glmer, type="pearson")^2))/df.residual(lh16glmer)

r.squaredGLMM(he15glmer)
r.squaredGLMM(lh14glmer)
r.squaredGLMM(lh15glmer)
r.squaredGLMM(lh16glmer)

roc.area(obs=HER15$for_bin,pred=fitted(he15glmer)) 
roc.area(obs=LHI14$for_bin,pred=fitted(lh14glmer)) 
roc.area(obs=LHI15$for_bin,pred=fitted(lh15glmer)) 
roc.area(obs=LHI16$for_bin,pred=fitted(lh16glmer)) 

#doesnt show anything like what I want.. try gam

he15gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(SQRTsmt, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=HER15,na.action=na.omit)

lh14gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(SQRTsmt, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=LHI14,na.action=na.omit)

lh15gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(SQRTsmt, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=LHI15,na.action=na.omit)

lh16gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(SQRTsmt, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=LHI16,na.action=na.omit)

sum((resid(he15gam, type="pearson")^2))/df.residual(he15gam)
sum((resid(lh14gam, type="pearson")^2))/df.residual(lh14gam)
sum((resid(lh15gam, type="pearson")^2))/df.residual(lh15gam)
sum((resid(lh16gam, type="pearson")^2))/df.residual(lh16gam)


roc.area(obs=HER15$for_bin,pred=fitted(he15gam)) 
roc.area(obs=LHI14$for_bin,pred=fitted(lh14gam)) 
roc.area(obs=LHI15$for_bin,pred=fitted(lh15gam)) 
roc.area(obs=LHI16$for_bin,pred=fitted(lh16gam)) 

# not going nowhere either..

glmer1<-glmer(for_bin~1+(1|trip_id),family="binomial", 
              data=HER15,na.action=na.omit)

HER15$glmer1_resid<-resid(glmer1)
p<-ggplot(data=HER15, aes(x=ekmU, y=glmer1_resid));p+geom_point()+stat_smooth(formula=y~s(x),method="gam")
summary(gam(glmer1_resid~s(ekmU, k=3), data=HER15))


summary(ghe15)
sum((resid(ghe15, type="pearson")^2))/df.residual(ghe15)

gam.check(ghe15)

roc.area(obs=HER15$for_bin,
         pred=fitted(ghe15)) 

