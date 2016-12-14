# 10/10/16 Lancaster, UK
# Model HMM data modelling foraging against non-foraging

rm(list=ls())
library(ggplot2)
library(reshape2)
library(lme4) 
library(verification)
library(MuMIn)
library(vegan)
library(ggrepel)
library(ncf)

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


setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("spreads/GPS_LT_141516_hmm_oceano_attribs.csv", h=T, strip.white=T)
# remove a few pesky LTs that have very little actual tracking data
# these are less than 24 hours long (< 144 rows)
#aggregate(HMMstates~trip_id, dat, FUN=function(x){length(x)})
dat<-dat[dat$trip_id != "18_02_16_21LW1"& dat$trip_id != "22_03_15_05RW2",]


d1<-melt(dat[,c(1:5, 26:44)], id.vars=c("Latitude", "Longitude", "Date", "Time", "TrackID"))

p<-ggplot(d1, aes(variable, value))
p+geom_boxplot()+facet_wrap(~variable, scale="free")

# investigate and remove outliers

dat[which(dat$AsstAG>20),] # !!!! WRAP in which() otherwise you get loads of NA rows!!!
#for moment we'll just turn to NA
dat[which(dat$AsstAG>20),]$AsstAG<-NA 

dat[which(dat$chlMH1>0.3),]
# these are reef-inflated CHL values, should probably leave rather than NA
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

#transformations change log to log2 and sqrt to smt/100 e.g. 100km unit change
dat$LOGCHL_M<-log2(dat$CHL_M)
dat$LOGCHL_A<-log2(dat$CHL_A)
dat$smt_100<-dat$smt/100

#collinearity
mcor<-cor(na.omit(dat[26:39]))
mcor 

library(car)
vif(lm(1:nrow(dat)~ekmU+modW+SST+ASST+
         sshHy+bty+smt_100+LOGCHL_A+yftA+yftJ+skjA+skjJ+
         betA+betJ,data=dat)) # all good

#start modelling
# get HMM into binary response
dat$for_bin<-0
dat[dat$HMMstates!=3,]$for_bin<-1

#setup NA removed col year datasets
# Very important to remove NA prior to modelling
d_mod<-dat[,c(6,20,16,23,30,31,36:42,44:46,50:52, 4, 5)]
HER15<-na.omit(d_mod[d_mod$Colony=="Heron",])
LHI14<-na.omit(d_mod[d_mod$Colony=="LHI" & d_mod$Year==2014,])
LHI15<-na.omit(d_mod[d_mod$Colony=="LHI"& d_mod$Year==2015,])
LHI16<-na.omit(d_mod[d_mod$Colony=="LHI"& d_mod$Year==2016,])

vif(lm(1:nrow(HER15)~ekmU+modW+SST+ASST+
         sshHy+bty+smt_100+LOGCHL_A+yftA+yftJ+skjA+skjJ+
         betA+betJ,data=HER15)) 

## RESAMPLE data to once every 3 points to reduce SPAC to acceptable level

HER15<-HER15[seq(1, nrow(HER15), 3),]
LHI14<-LHI14[seq(1, nrow(LHI14), 3),]
LHI15<-LHI15[seq(1, nrow(LHI15), 3),]
# LHI16 doesnt have tuna data so dropping

#removeunused factors
HER15$trip_id<-factor(HER15$trip_id)
LHI14$trip_id<-factor(LHI14$trip_id)
LHI15$trip_id<-factor(LHI15$trip_id)

# do z transform if needed due to convergence error - often corr. related
#library(vegan)
#HER15<-data.frame(HER15[,c(18,1,2,3,19, 20)],decostand(HER15[,4:17],method="standardize"))
#LHI14<-data.frame(LHI14[,c(12,1,2,3)],decostand(LHI14[,4:11],method="standardize"))
#LHI15<-data.frame(LHI15[,c(12,1,2,3)],decostand(LHI15[,4:11],method="standardize"))

#**** ^^^%%%% Heron Island %%%^^^ *****#
# corr
pairs(HER15[,4:17 ], upper.panel = panel.smooth,lower.panel=panel.cor)
# outliers, in tuna data especially
cor(HER15[,4:17 ])

enviro_std<-decostand(HER15[,4:17 ], method="standardize")
# takes all varibs but eke also juv and adu tnua forms

#!! then do pca (just a scaled RDA) !!#
enviro_rda<-rda(enviro_std, scale=T)
summary(enviro_rda, display=NULL)
screeplot(enviro_rda) # badly scaled
#full summary
#summary(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1)) 
# i've put scaling to 1 for the sites to fit better on the plot
# Now make some plots
enviro.species.scores<-as.data.frame(scores(enviro_rda, display='species'))
enviro.species.scores$Predictors<-colnames(enviro_std)
#enviro.species.scores$Pred_codes<-codez
head(enviro.species.scores)

g<- ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_point(data=enviroPCA, aes(y=PC2, x=PC1, shape=treatshape, fill=treatfill),size=3)+scale_shape_identity()+scale_fill_identity()+
  geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=PC2, xend=PC1), arrow=arrow(length=unit(0.3,'lines')), colour='red')+theme_classic() 
g<-g+geom_text_repel(data=enviro.species.scores, aes(y=PC2, x=PC1, label=Predictors), segment.size=0, colour='red')

eig<-eigenvals(enviro_rda)
g<- g+scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))

g
###!! Save RDA !!### 

# Have a look at data
d1<-melt(HER15, id.vars=c("for_bin","trip_id", "Colony", "Year"))
g1<-ggplot(data=d1, aes(x=value))
g1+geom_histogram()+facet_wrap(~variable, scale="free")
# look to remove outliers
#HER15<-HER15[-which(HER15$skjA>0.1),]
d1<-melt(HER15, id.vars=c("for_bin","trip_id", "Colony", "Year"))
g1<-ggplot(data=d1, aes(x=value))
g1+geom_histogram()+facet_wrap(~variable, scale="free")

g1<-ggplot(data=d1, aes(y=value, x=factor(for_bin)))
g1+geom_boxplot()+facet_wrap(~variable, scale="free")

g1+geom_jitter(height=0.1)+geom_smooth(method="glm")+facet_wrap(~variable, scale="free")

# selecting which tuna variables to use
AIC(glmer(for_bin~skjA + (1|trip_id), family="binomial", data=HER15),
    glmer(for_bin~skjJ + (1|trip_id), family="binomial", data=HER15),
    glmer(for_bin~yftA + (1|trip_id), family="binomial", data=HER15),
    glmer(for_bin~yftJ + (1|trip_id), family="binomial", data=HER15),
    glmer(for_bin~betA + (1|trip_id), family="binomial", data=HER15),
    glmer(for_bin~betJ + (1|trip_id), family="binomial", data=HER15))
 
# yft adu and bet juv   
pairs(HER15[,c(4:7, 12, 15,17) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# not bad, try without modW 

qplot(data=HER15, x=yftA, bins=50)+facet_grid(for_bin~.)

# Easier to rescale variables after analysis
#bet_juvenil_potential_biomass
#BIGEYE weekly biomass distribution
#HER15$yftA<-(HER15$yftA*1000000)/1000 
#HER15$betJ<-(HER15$betJ*1000000)/1000 

# Choose which variables are better suited to 2nd degree polynomial
g1<-ggplot(data=d1, aes(y=for_bin, x=value))
g1+geom_jitter(height=0.1)+
  geom_smooth(method="glm")+
  geom_smooth(method="glm", formula=y~poly(x, 2), colour=2)+
  facet_wrap(~variable, scale="free")
#polys for ekm and ssh

he15glmer<-glmer(for_bin~poly(ekmU,2)+modW+ASST+poly(sshHy,2)+
            smt_100+ yftA + betJ+(1|trip_id), family="binomial", 
            data=HER15)
# throws warning this is due to corr variables and yft being crap
# can get round by rescaling, in any case YFT could maybe go
summary(he15glmer)
sum((resid(he15glmer, type="pearson")^2))/df.residual(he15glmer)

# test if yft can go
he15glmer2<-glmer(for_bin~poly(ekmU,2)+modW+ASST+poly(sshHy,2)+
                   smt_100+ betJ+(1|trip_id), family="binomial", 
                 data=HER15)

summary(he15glmer2) # modW looking reversed. due to corr remove

he15glmer3<-glmer(for_bin~poly(ekmU,2)+ASST+poly(sshHy,2)+
                    smt_100+ betJ+(1|trip_id), family="binomial", 
                  data=HER15)

anova(he15glmer,he15glmer2, he15glmer3) # drop it and get rid of wnd
confint(he15glmer3, method='boot', nsim=99) # takes ages
drop1(he15glmer2, test="Chisq") # no others

# Model details

sum((resid(he15glmer3, type="pearson")^2))/df.residual(he15glmer3)
#1.023228
summary(he15glmer3)
roc.area(obs=HER15$for_bin,pred=fitted(he15glmer3)) #0.8
r.squaredGLMM(he15glmer3)
#R2m       R2c 
#0.3609434 0.4977527 


ms1<-model.sel(glmer(for_bin~smt_100+ betJ+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~poly(ekmU,2)+ betJ+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~modW+betJ+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~ASST+betJ+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~poly(sshHy,2)+betJ+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~betJ+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~smt_100+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~poly(ekmU,2)+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~poly(sshHy,2)+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~ASST+(1|trip_id), family="binomial", data=HER15),
          glmer(for_bin~modW+(1|trip_id), family="binomial", data=HER15))

# ok interesting for sure rank via AIC or BIC, doesnt make much diff
importance(ms1) # very important number of times term is used are equal

# Either go with the above AIC comp or just dont do importance all together
#re. http://blog.minitab.com/blog/adventures-in-statistics/how-to-identify-the-most-important-predictor-variables-in-regression-models

# check spatial autocorr of final model, should be fine due to 1:3 point resample
corm3 <- spline.correlog(HER15$Longitude,
         HER15$Latitude, residuals(he15glmer2, type="pearson"), 
         na.rm=T, latlon=T,resamp=1)

#**** ^^^%%%% %%%^^^ *****#

#**** ^^^%%%% LHI 2014 %%%^^^ *****#
# corr
pairs(LHI14[,5:18 ], upper.panel = panel.smooth,lower.panel=panel.cor)
# outliers, in tuna data especially
cor(LHI14[,5:18 ])

enviro_std<-decostand(LHI14[,5:18 ], method="standardize")
# takes all varibs but eke also juv and adu tnua forms

#!! then do pca (just a scaled RDA) !!#
enviro_rda<-rda(enviro_std, scale=T)
summary(enviro_rda, display=NULL)
screeplot(enviro_rda) # badly scaled
#full summary
#summary(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1)) 
# i've put scaling to 1 for the sites to fit better on the plot
# Now make some plots
enviro.species.scores<-as.data.frame(scores(enviro_rda, display='species'))
enviro.species.scores$Predictors<-colnames(enviro_std)
#enviro.species.scores$Pred_codes<-codez
head(enviro.species.scores)

g<- ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_point(data=enviroPCA, aes(y=PC2, x=PC1, shape=treatshape, fill=treatfill),size=3)+scale_shape_identity()+scale_fill_identity()+
  geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=PC2, xend=PC1), arrow=arrow(length=unit(0.3,'lines')), colour='red')+theme_classic() 
g<-g+geom_text_repel(data=enviro.species.scores, aes(y=PC2, x=PC1, label=Predictors), segment.size=0, colour='red')

eig<-eigenvals(enviro_rda)
g<- g+scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))

g
###!! Save RDA !!### 

# Have a look at data
d1<-melt(LHI14, id.vars=c("for_bin","trip_id", "Colony", "Year"))
g1<-ggplot(data=d1, aes(x=value))
g1+geom_histogram()+facet_wrap(~variable, scale="free")
# look to remove outliers
#LHI14<-LHI14[-which(LHI14$skjA>0.1),]
d1<-melt(LHI14, id.vars=c("for_bin","trip_id", "Colony", "Year", "Month"))
g1<-ggplot(data=d1, aes(x=value))
g1+geom_histogram()+facet_wrap(~variable, scale="free")

g1<-ggplot(data=d1, aes(y=value, x=factor(for_bin)))
g1+geom_boxplot()+facet_wrap(~variable, scale="free")

g1+geom_jitter(height=0.1)+geom_smooth(method="glm")+facet_wrap(~variable, scale="free")

# selecting which tuna variables to use
AIC(glmer(for_bin~skjA + (1|trip_id), family="binomial", data=LHI14),
    glmer(for_bin~skjJ + (1|trip_id), family="binomial", data=LHI14),
    glmer(for_bin~yftA + (1|trip_id), family="binomial", data=LHI14),
    glmer(for_bin~yftJ + (1|trip_id), family="binomial", data=LHI14),
    glmer(for_bin~betA + (1|trip_id), family="binomial", data=LHI14),
    glmer(for_bin~betJ + (1|trip_id), family="binomial", data=LHI14))

# yft juv and bet juv   + skj adu
pairs(LHI14[,c(5:7,9, 13, 14,16,17,18)], upper.panel = panel.smooth,lower.panel=panel.cor)
# not bad, try without modW 


qplot(data=LHI14, x=LOGCHL_A, bins=50)+facet_grid(for_bin~.)
# remove CHL outliers
LHI14<-LHI14[which(LHI14$LOGCHL_A< -2),]

# Easier to rescale variables after analysis
#bet_juvenil_potential_biomass
#BIGEYE weekly biomass distribution
#LHI14$yftA<-(LHI14$yftA*1000000)/1000 
#LHI14$betJ<-(LHI14$betJ*1000000)/1000 

# Choose which variables are better suited to 2nd degree polynomial
g1<-ggplot(data=d1, aes(y=for_bin, x=value))
g1+geom_jitter(height=0.1)+
  geom_smooth(method="glm")+
  geom_smooth(method="glm", formula=y~poly(x, 2), colour=2)+
  facet_wrap(~variable, scale="free")
#polys for ekm and ssh

lhi14glmer<-glmer(for_bin~poly(ekmU,2)+modW+poly(ASST,2)+poly(sshHy,2)+
                   smt_100+ yftJ + LOGCHL_A+ (1|trip_id), family="binomial", 
                 data=LHI14)
# throws warning this is due to corr variables and smt being crap
# can get round by rescaling, in any case modW could maybe go
summary(lhi14glmer)
sum((resid(lhi14glmer, type="pearson")^2))/df.residual(lhi14glmer)

# We remove modW and refit - still warning 
lhi14glmer2<-glmer(for_bin~poly(ekmU,2)+poly(ASST, 2)+poly(sshHy,2)+
                    smt_100+ yftJ + LOGCHL_A+(1|trip_id), family="binomial", 
                  data=LHI14)

# We remove smt and refit - no warning
lhi14glmer3<-glmer(for_bin~poly(ekmU,2)+poly(ASST,2)+poly(sshHy,2)+
                    yftJ + LOGCHL_A+(1|trip_id), family="binomial", 
                   data=LHI14)

anova(lhi14glmer,lhi14glmer2, lhi14glmer3) # good to drop modW and smt
summary(lhi14glmer2) # some other variables look iffy now, test with confint

confint(lhi14glmer2, method='boot', nsim=99) # takes ages remove smt and ssh poly
confint(lhi14glmer2)# no warning
drop1(lhi14glmer3, test="Chisq") # Asst and ssh's poly?

lhi14glmer3a<-glmer(for_bin~poly(ekmU,2)+ASST+poly(sshHy,2)+
                     yftJ + LOGCHL_A+(1|trip_id), family="binomial", 
                   data=LHI14)
lhi14glmer3b<-glmer(for_bin~poly(ekmU,2)+poly(ASST,2)+sshHy+
                      yftJ + LOGCHL_A+(1|trip_id), family="binomial", 
                    data=LHI14)
lhi14glmer3c<-glmer(for_bin~poly(ekmU,2)+ASST+sshHy+
                      yftJ + LOGCHL_A+(1|trip_id), family="binomial", 
                    data=LHI14)

anova(lhi14glmer3,lhi14glmer3a, lhi14glmer3b, lhi14glmer3c) 
# poly terms not significant
lhi14glmer3<-lhi14glmer3c
# Model details

sum((resid(lhi14glmer3, type="pearson")^2))/df.residual(lhi14glmer3)
# 0.9804307
summary(lhi14glmer3)
roc.area(obs=LHI14$for_bin,pred=fitted(lhi14glmer3)) #0.73
r.squaredGLMM(lhi14glmer3)
#R2m       R2c 
#0.1492239 0.2851388   

# check spatial autocorr of final model, should be fine due to 1:3 point resample
corm3 <- spline.correlog(LHI14$Longitude,
                         LHI14$Latitude, residuals(lhi14glmer3, type="pearson"), 
                         na.rm=T, latlon=T,resamp=1)

#**** ^^^%%%% %%%^^^ *****#

#**** ^^^%%%% LHI 2014 %%%^^^ *****#
# corr
pairs(LHI15[,5:18 ], upper.panel = panel.smooth,lower.panel=panel.cor)
# outliers, in tuna data especially
cor(LHI15[,5:18 ])

enviro_std<-decostand(LHI15[,5:18 ], method="standardize")
# takes all varibs but eke also juv and adu tnua forms

#!! then do pca (just a scaled RDA) !!#
enviro_rda<-rda(enviro_std, scale=T)
summary(enviro_rda, display=NULL)
screeplot(enviro_rda) # badly scaled
#full summary
#summary(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1)) 
# i've put scaling to 1 for the sites to fit better on the plot
# Now make some plots
enviro.species.scores<-as.data.frame(scores(enviro_rda, display='species'))
enviro.species.scores$Predictors<-colnames(enviro_std)
#enviro.species.scores$Pred_codes<-codez
head(enviro.species.scores)

g<- ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_point(data=enviroPCA, aes(y=PC2, x=PC1, shape=treatshape, fill=treatfill),size=3)+scale_shape_identity()+scale_fill_identity()+
  geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=PC2, xend=PC1), arrow=arrow(length=unit(0.3,'lines')), colour='red')+theme_classic() 
g<-g+geom_text_repel(data=enviro.species.scores, aes(y=PC2, x=PC1, label=Predictors), segment.size=0, colour='red')

eig<-eigenvals(enviro_rda)
g<- g+scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))

g
###!! Save RDA !!### 

# Have a look at data
d1<-melt(LHI15, id.vars=c("for_bin","trip_id", "Colony", "Year"))
g1<-ggplot(data=d1, aes(x=value))
g1+geom_histogram()+facet_wrap(~variable, scale="free")
# look to remove outliers
#LHI15<-LHI15[-which(LHI15$skjA>0.1),]
d1<-melt(LHI15, id.vars=c("for_bin","trip_id", "Colony", "Year", "Month"))
g1<-ggplot(data=d1, aes(x=value))
g1+geom_histogram()+facet_wrap(~variable, scale="free")

g1<-ggplot(data=d1, aes(y=value, x=factor(for_bin)))
g1+geom_boxplot()+facet_wrap(~variable, scale="free")

# selecting which tuna variables to use
AIC(glmer(for_bin~skjA + (1|trip_id), family="binomial", data=LHI15),
    glmer(for_bin~skjJ + (1|trip_id), family="binomial", data=LHI15),
    glmer(for_bin~yftA + (1|trip_id), family="binomial", data=LHI15),
    glmer(for_bin~yftJ + (1|trip_id), family="binomial", data=LHI15),
    glmer(for_bin~betA + (1|trip_id), family="binomial", data=LHI15),
    glmer(for_bin~betJ + (1|trip_id), family="binomial", data=LHI15))

# skj adu or juv   + yft juv
pairs(LHI15[,c(5:7,10,11, 14,16,17,18)], upper.panel = panel.smooth,lower.panel=panel.cor)
# not bad, try without modW 


qplot(data=LHI15, x=ASST, bins=50)+facet_grid(for_bin~.)

# Easier to rescale variables after analysis
#bet_juvenil_potential_biomass
#BIGEYE weekly biomass distribution
#LHI15$yftA<-(LHI15$yftA*1000000)/1000 
#LHI15$betJ<-(LHI15$betJ*1000000)/1000 

# Choose which variables are better suited to 2nd degree polynomial
g1<-ggplot(data=d1, aes(y=for_bin, x=value))
g1+geom_jitter(height=0.1)+
  geom_smooth(method="glm")+
  geom_smooth(method="glm", formula=y~poly(x, 2), colour=2)+
  facet_wrap(~variable, scale="free")
#polys for ekm and ASST, can't have ssh

LHI15glmer<-glmer(for_bin~poly(ekmU,2)+modW+poly(ASST,2)+
                    smt_100+ skjA+  (1|trip_id), family="binomial", 
                  data=LHI15)

summary(LHI15glmer) # modW looks inflated, probs due to corr
sum((resid(LHI15glmer, type="pearson")^2))/df.residual(LHI15glmer)

# We remove modW and refit 
LHI15glmer2<-glmer(for_bin~poly(ekmU,2)+poly(ASST,2)+
                     smt_100+ skjA+  (1|trip_id), family="binomial", 
                   data=LHI15)


anova(LHI15glmer,LHI15glmer2) # higher AIC in glmer2 but need to drop wnd
summary(LHI15glmer2) # good, smt looks less sig tho

confint(LHI15glmer2, method='boot', nsim=99) # takes ages remove smt and ssh poly
confint(LHI15glmer2)# no warning
drop1(LHI15glmer2, test="Chisq") # sqrt?

#remove smt

LHI15glmer2<-glmer(for_bin~poly(ekmU,2)+poly(ASST,2)+
                     skjA+  (1|trip_id), family="binomial", 
                   data=LHI15)


# Model details

sum((resid(LHI15glmer2, type="pearson")^2))/df.residual(LHI15glmer2)
#1.068224
summary(LHI15glmer2)
roc.area(obs=LHI15$for_bin,pred=fitted(LHI15glmer2)) #0.73
r.squaredGLMM(LHI15glmer2)
#R2m       R2c 
#0.2965124 0.5212758  

# check spatial autocorr of final model, should be fine due to 1:3 point resample
corm3 <- spline.correlog(LHI15$Longitude,
                         LHI15$Latitude, residuals(LHI15glmer2, type="pearson"), 
                         na.rm=T, latlon=T,resamp=1)





#### old

lh14glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
                   smt_100+LOGCHL_A + (1|trip_id), family="binomial", 
                 data=LHI14,na.action=na.omit)
summary(lh14glmer)
sum((resid(lh14glmer, type="pearson")^2))/df.residual(lh14glmer)

lh15glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
                   smt_100+LOGCHL_A + (1|trip_id), family="binomial", 
                 data=LHI15,na.action=na.omit)
summary(lh15glmer)
sum((resid(lh15glmer, type="pearson")^2))/df.residual(lh15glmer)

lh16glmer<-glmer(for_bin~ekmU+modW+SST+ASST+sshHy+bty+
                   smt_100+LOGCHL_A + (1|trip_id), family="binomial", 
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
               s(smt_100, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=HER15,na.action=na.omit)

lh14gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(smt_100, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=LHI14,na.action=na.omit)

lh15gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(smt_100, k=5)+s(LOGCHL_A, k=5), family="binomial", 
             data=LHI15,na.action=na.omit)

lh16gam<-gam(for_bin~s(ekmU, k=5)+s(modW, k=5)+s(SST, k=5)+s(ASST, k=5)+
               s(sshHy, k=5)+s(bty, k=5)+
               s(smt_100, k=5)+s(LOGCHL_A, k=5), family="binomial", 
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

### ~~~ VISUALISING DAY NIGHT stuff ~~~ ###

out<-NULL
for(i in unique(dat$trip_id))
{
  d0<-with(dat[dat$trip_id==i,], tapply(HMMstates, list(Hr, HMMstates), FUN = length))
  d1<-data.frame(trip_id=i, Hr=dimnames(d0)[[1]], d0, Colony=unique(dat[dat$trip_id==i,]$Colony), Year=unique(dat[dat$trip_id==i,]$Year))
  out<-rbind(out, d1)
  print(i)
}
out[is.na(out)] <- 0
out$hmmSUM<-rowSums(out[,3:5], na.rm=T)
out$propsit<-round(out$X1/out$hmmSUM*100)
out$propfrg<-round(out$X2/out$hmmSUM*100)
out$propfly<-round(out$X3/out$hmmSUM*100)

#for plotting

d1<-melt(out[,c(1:2, 6,7,9:11)], id.vars=c("trip_id", "Hr", "Colony", "Year"))

g<-ggplot(data=d1[d1$Colony=="Heron",], aes(y=value, x=factor(Hr), fill=variable))
g+geom_boxplot()+xlim(paste(seq(0,23)))+theme_classic()

g<-ggplot(data=d1[d1$Colony=="LHI"&d1$Year==2014,], aes(y=value, x=factor(Hr), fill=variable))
g+geom_boxplot()+xlim(paste(seq(0,23)))+theme_classic()

g<-ggplot(data=d1[d1$Colony=="LHI"&d1$Year==2015,], aes(y=value, x=factor(Hr), fill=variable))
g+geom_boxplot()+xlim(paste(seq(0,23)))+theme_classic()

g<-ggplot(data=d1[d1$Colony=="LHI"&d1$Year==2016,], aes(y=value, x=factor(Hr), fill=variable))
g+geom_boxplot()+xlim(paste(seq(0,23)))+theme_classic()

d1$c1<-paste(d1$Colony, d1$Year)

g<-ggplot(data=d1, aes(y=value, x=factor(Hr), fill=variable))
g+geom_boxplot(outlier.size=0)+xlim(paste(seq(0,23)))+theme_classic()+facet_wrap(~c1)

g<-ggplot(data=d1, aes(y=value, x=factor(Hr), fill=variable))
g+geom_boxplot(outlier.size=0)+xlim(paste(seq(0,23)))+theme_classic()+facet_wrap(~Colony)

g<-ggplot(data=d1,aes(y=value, x=Hr, group=variable))
g+geom_jitter(aes(colour=variable), width= 0.2, size=1)+
  stat_summary(fun.y=mean, geom="line", size=2,  aes(color=variable))+
  xlim(paste(seq(0,23)))+theme_classic()+facet_wrap(~c1)

# hmm just stick to summary tbles

#experimental density plot

d2<-NULL
for(i in 1: nrow(d1)){
  d20<-data.frame(Colony=rep(d1[i,]$Colony, d1[i,]$value),
                  Year=rep(d1[i,]$Year, d1[i,]$value),
                  Hr=rep(d1[i,]$Hr, d1[i,]$value), 
                  Behaviour=rep(d1[i,]$variable, d1[i,]$value))
  d2<-rbind(d2, d20)
} # takes a while!
ggplot(d2, aes(x=as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
  geom_histogram(bins=24)

ggplot(d2[d2$Colony=="Heron",], aes(as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
  geom_density(alpha = 0.1)+xlim(paste(seq(0,23)))+theme_classic()
# note as.numeric(as.character(Hr) - has to be or factors screw it up
# might be able to improve 
ggplot(d2[d2$Colony=="LHI" & d2$Year==2014,], aes(as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
  geom_density(alpha = 0.1)+xlim(paste(seq(0,23)))+theme_classic()

ggplot(d2[d2$Colony=="LHI" & d2$Year==2015,], aes(as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
  geom_density(alpha = 0.1)+xlim(paste(seq(0,23)))+theme_classic()

ggplot(d2[d2$Colony=="LHI" & d2$Year==2016,], aes(as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
  geom_density(alpha = 0.1)+xlim(paste(seq(0,23)))+theme_classic()

d2$c1<-paste(d2$Colony, d2$Year)

ggplot(d2, aes(as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
  geom_density(alpha = 0.1)+xlim(paste(seq(0,23)))+theme_classic()+facet_wrap(~c1)

#Attempt to out night day stuff in there
#ggplot(d2, aes(as.numeric(as.character(Hr)), fill = Behaviour, colour = Behaviour)) +
#  geom_rect(aes(xmin = 0, xmax = 6, ymin=0, ymax=0.2, colour="grey" ))+
#  geom_rect(aes(xmin = 18, xmax =23 , ymin=0, ymax=0.2, colour="grey"  ))+
#  geom_density(alpha = 0.1)+xlim(paste(seq(0,23)))+theme_classic()+facet_wrap(~Colony)

# Now day night table

out<-NULL
for(i in unique(dat$trip_id))
{
  d0<-with(dat[dat$trip_id==i,], tapply(HMMstates, list(DayNight, HMMstates), FUN = length))
  d1<-data.frame(trip_id=i, DN=dimnames(d0)[[1]], d0, Colony=unique(dat[dat$trip_id==i,]$Colony), Year=unique(dat[dat$trip_id==i,]$Year))
  out<-rbind(out, d1)
  print(i)
}
out[is.na(out)] <- 0
out$hmmSUM<-rowSums(out[,3:5], na.rm=T)
out$propsit<-round(out$X1/out$hmmSUM*100)
out$propfrg<-round(out$X2/out$hmmSUM*100)
out$propfly<-round(out$X3/out$hmmSUM*100)

d1<-melt(out[,c(1:2, 6,7,9:11)], id.vars=c("trip_id", "DN", "Colony", "Year"))


out_DN<-rbind(data.frame(ID="Heron", with(d1[d1$Colony=="Heron",],
                                          tapply(value, list(DN, variable), FUN = mean))),
              data.frame(ID="Heron",with(d1[d1$Colony=="Heron",],
                                         tapply(value, list(DN, variable), FUN = sd))),
              data.frame(ID="LHI15", with(d1[d1$Colony=="LHI"& out$Year==2015,],
                                          tapply(value, list(DN, variable), FUN = mean))),
              data.frame(ID="LHI15", with(d1[d1$Colony=="LHI"& out$Year==2015,],
                                          tapply(value, list(DN, variable), FUN = sd))),
              data.frame(ID="LHI14", with(d1[d1$Year==2014,],
                                          tapply(value, list(DN, variable), FUN = mean))),
              data.frame(ID="LHI14", with(d1[d1$Year==2014,],
                                          tapply(value, list(DN, variable), FUN = sd))),
              data.frame(ID="LHI16", with(d1[d1$Year==2016,],
                                          tapply(value, list(DN, variable), FUN = mean))),
              data.frame(ID="LHI16", with(d1[d1$Year==2016,],
                                          tapply(value, list(DN, variable), FUN = sd))))

write.csv(out_DN, "paper_results/DayNight_HMM_summary.csv" , quote=F, row.names=T)

