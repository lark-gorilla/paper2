# 19/11/16 LBishops Castle, UK
# Model climatology data comparing foraging areas against background 
# Final paper version including method to deal with SPAC

rm(list=ls())
library(ggplot2)
library(ggrepel)
library(reshape2)
library(sp)
library(rgeos)
library(car)
library(vegan)
library(verification)
library(spdep)
library(ncf)
library(pscl)
library(survey)
library(caret)
#install_github("clbustos/dominanceAnalysis")
library(dominanceanalysis)
source("~/grive/phd/scripts/paper2/hackedDominance.R")

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
#dat[is.na(dat$chl),]$wnd<-NA
# and some weird sst
dat[which(dat$sst<0),]$sst<-NA
#really slow winds near png<-NA
#dat[which(dat$wnd<6),]$wnd<-NA
#weird ekem point
#dat[which(dat$ekm< -0.00001),]$ekm<-NA
#remove skj_juv outliers
#dat[which(dat$skj_juv_03_ave>28),]$skj_juv_03_ave<-NA

dat$PA<-0
dat[dat$Count>0,]$PA<-1 # making presence absence variable

# remove some outlier presence clusters that cause modeeling issues later
dat<-dat[-which(dat$dset=="Heron"& dat$PA==1 & dat$Longitude<153.04 &
                 dat$Latitude< -21.44),]

dat<-dat[-which(dat$dset=="Heron"& dat$PA==1 & dat$Longitude<152.5 &
                 dat$Latitude< -18.74),]

dat<-dat[-which(dat$dset=="Heron" & dat$PA==1 & dat$Longitude<151.4),]

dat<-dat[-which(dat$dset=="Heron" & dat$PA==1 & dat$Longitude<153.953 &
                  dat$Latitude< -25.358),]

dat<-dat[-which(dat$dset=="Heron" & dat$PA==1 & dat$chl>1000),]

## NOW we do naomit! 
## This is integral to modelling with the IT approach and others, because
## when subset models are created, they have different numbers of NAs
## and therefore data: they are not comparable

nrow(dat)
nrow(na.omit(dat))
dat<-na.omit(dat)

#write.csv(dat, "test4.csv", quote=F, row.names=F) # write out na omitted dat

# seamounts and chl still need some attention
qplot(chl, data=dat);qplot(log(chl), data=dat)
qplot(smt, data=dat);qplot(sqrt(smt), data=dat)

# Transform variables AFTER naomit!!
dat$chl_log<-log2(dat$chl)
dat$smt_100-dat$smt/100
# changes BET from g/m2 to kg/km2 (bigger units in model)
#dat$bet_adu_03_ave<-(dat$bet_adu_03_ave*1000000)/1000 

qplot(data=dat[dat$dset=="Heron",], x=wnd, bins=50)+facet_grid(PA~.)

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
heronPA<-heronA[sample(1:nrow(heronA), 3200),]
dat_heron<-rbind(heronP, heronPA)

# Below code removed duplicate points where both presence and absence occur
# together, however caused perfect seperation in binomial glm
#nrow(dat_heron)   
#sp1<-SpatialPointsDataFrame(dat_heron[,1:2], data=dat_heron)
#sp2<-remove.duplicates(sp1, remove.second = TRUE, memcmp = TRUE) # removes duplicate points, remove.second makes sure presences remain and absences are removed
#dat_heron<-sp2@data
#nrow(dat_heron) 

# finding distance to nearest point after na.omit, this will identify isolated points
# proj dat into m units for gDistance
sp1<-SpatialPointsDataFrame(dat_heron[,1:2], data=dat_heron, proj4string=CRS("+proj=longlat + datum=wgs84"))
sp1proj <- spTransform(sp1, CRS=CRS("+proj=laea +lon_0=155 +lat_0=-21"))

d <- gDistance(sp1proj, byid=T)

min.d <- apply(d, 1, function(x) min(x[x>1000])) # pulls nearest point, which isnt itself or within 1000m projection error range
dat_heron$min.d<-min.d
hist(dat_heron$min.d)

plot(Latitude~Longitude, dat_heron)
points(Latitude~Longitude, dat_heron[dat_heron$min.d>50000,], col=2)

dat_heron<-dat_heron[dat_heron$min.d<50000,]



# now have a look at collinearity
vif(lm(1:nrow(dat_heron)~sst+shd+ekm+chl_log+wnd+tmc+smt_100+bty+bet_adu_03_ave+
         bet_juv_03_ave+bet_tot_03_ave+yft_adu_03_ave+
         yft_juv_03_ave+yft_tot_03_ave+skj_adu_03_ave+
         skj_juv_03_ave+skj_tot_03_ave, data=dat_heron))  

vif(lm(1:nrow(dat_heron)~shd+ekm+chl_log+wnd+tmc+smt_100+bty+bet_adu_03_ave+
         bet_juv_03_ave+yft_adu_03_ave+
         skj_adu_03_ave,
       data=dat_heron))

# ok looks good, all values < 10, hmm skjadu and yftadu are corr according to 
cor(dat_heron[,-(1:4)]) # NOTE: 0.7 seems far to 

#pairs(dat_heron[,-(1:2)], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_heron[,5:23], upper.panel = panel.smooth,lower.panel=panel.cor)
# ok so actually for Heron I include both BET then
# adult yft and juv skj
# final one for model
pairs(dat_heron[,c(5,6,9,10,12,13,14,16,17, 19,20, 23, 24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# looking at it in full, i would be tempted to remove shd, tmc and yellowfin adult

# Some serious collinearity need to investigate using PCA
#then z transform to standardize
enviro_std<-decostand(dat_heron[,c(5,6,9,10,12,13,14,16,17, 19,20, 23, 24)], method="standardize")
# takes all varibs but eke also juv and adu tnua forms

# then do pca (just a scaled RDA)
enviro_rda<-rda(enviro_std, scale=T)
summary(enviro_rda, display=NULL)
screeplot(enviro_rda) # badly scaled
#full summary
summary(enviro_rda)

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

# Basically for Heron we have ~50% of the variation explained by PC1 
# Which is effectively a gradient of Latitude with increasing SST and SKJ
# and YFT numbers (Coral Sea) at one end and High Chlorophyll, SHD 
# and deep thermocline at the other end (Tasman sea)

summary(enviro_rda)$species

glm_pca<-data.frame(PA=dat_heron$PA, enviro.sites.scores)
##names(glm_pca)[names(glm_pca)=="PC1"]<-"inv_latitude"

#try PC regression
m2<-glm(PA~PC1+PC2+PC3+PC4, data=glm_pca, family="binomial")

resglm<-residuals(m2, type="pearson")
summary(m2)
print(sum((resid(m2, type="pearson")^2))/df.residual(m2))
print(roc.area(dat_heron$PA, fitted(m2))$A)


# check final correltaions and VIFs
pairs(dat_heron[,c(5,9,12,13 ,23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)

pairs(dat_heron[,c(9,10,12,13 ,23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# 

vif(lm(1:nrow(dat_heron)~bet_adu_03_ave+
        bty+smt_100+wnd+sst+chl_log,
       data=dat_heron))

# Now we see which variables are better suited to a second deg polynomial

d1<-melt(dat_heron, id.vars=c("PA", "Count", "dset"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=PA, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")

# We would do below for each variable individually to see the suitability of a poly

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
# not using ekm as its crap anyway

# To account for spatial autocorrelation we calculate autocovariate terms
# using the residuls of the fitted model and a neighbourhood distance to 
# link spatially correlated variables, this term (RAC) is then included in
# the final model to soak up the residual deviance 

# Choose either sst, chla, skj_adu as tropical tuna metric

# see which tuna metric is best predictor

AIC(glm(PA~bet_adu_03_ave,data=dat_heron, family="binomial"),
    glm(PA~bet_juv_03_ave,data=dat_heron, family="binomial"),
    glm(PA~skj_adu_03_ave,data=dat_heron, family="binomial"),
    glm(PA~skj_juv_03_ave,data=dat_heron, family="binomial"),
    glm(PA~yft_adu_03_ave,data=dat_heron, family="binomial"),
    glm(PA~yft_juv_03_ave,data=dat_heron, family="binomial"))

# bet_juv, skj_juv, bet_adu, skj_adu    
  
pairs(dat_heron[,c(12,13,14,16,17, 20 ,23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)

pairs(dat_heron[,c(9,12,13,16 ,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_heron[,c(5,9,12,13 ,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_heron[,c(9,12,13, 23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)

#outlers again
qplot(data=dat_heron, x=chl_log, bins=50)+facet_grid(PA~.)

dat_heron<-dat_heron[dat_heron$bet_adu_03_ave<0.012,]
dat_heron<-dat_heron[dat_heron$chl_log<7,]

m2a<-glm(PA~wnd+chl_log+smt_100+bty+tmc,
  data=dat_heron, family="binomial")
m2b<-glm(PA~wnd+sst+smt_100+bty,
  data=dat_heron, family="binomial")
m2c<-glm(PA~skj_adu_03_ave+bet_juv_03_ave+
           smt_100+bty+bet_adu_03_ave,
         data=dat_heron, family="binomial")
m2d<-glm(PA~wnd+shd+smt_100+bty+tmc,
         data=dat_heron, family="binomial")

AIC(m2a, m2b, m2c, m2d)


anova(m2a, m2b, m2c);AIC(m2a, m2b, m2c)
pR2(m2a)[4];pR2(m2b)[4];pR2(m2c)[4]

m2<-glm(PA~chl_log+smt_100+yft_juv_03_ave+
          bty+bet_juv_03_ave,
         data=dat_heron, family="binomial");resglm<-residuals(m2, type="pearson")

resglm<-residuals(m2, type="pearson")
sp2<-SpatialPointsDataFrame(dat_heron[,1:2], data=data.frame(resglm), 
    proj4string=CRS("+proj=longlat + datum=wgs84"))
bubble(sp2, zcol='resglm')
#bubble(sp2, zcol='resglm', do.sqrt=F)


#write.csv(cbind(dat_heron, resglm), "test5.csv", quote=FALSE)


#corglm <- correlog(dat_heron$Longitude, dat_heron$Latitude, residuals(m2, type="pearson"), na.rm=T,
#                   latlon=T, increment=25,resamp=1)
RAC<-autocov_dist(resglm, cbind(dat_heron[,1], dat_heron[,2]),
                  nbs = 50, type = "one", zero.policy = T,
                  style = "B", longlat=TRUE)
# Model outputs are sensitive to neighborhood distance.
# here we use 50 kms as happy medium between different
# variables' resolutions. Remember points are at ~ 10km
# as are chl and bty smt tuna are 100km resampled at 25.
# 50km neighbourhood gives better SPAC reduction than 75

dat_heron$RAC<-RAC

m3<-glm(PA~chl_log+smt_100+yft_juv_03_ave+
          bty+bet_juv_03_ave+RAC,
        data=dat_heron, family="binomial")

plot(m3)
print(sum((resid(m3, type="pearson")^2))/df.residual(m3))

# Huge resids = overdisp.
# refit on reidual removed dataset

dat_heron<-dat_heron[-which(resid(m3, type="pearson")< -10),]

m2<-glm(PA~chl_log+smt_100+yft_juv_03_ave+
          bty+bet_juv_03_ave,
        data=dat_heron, family="binomial");resglm<-residuals(m2, type="pearson")

RAC<-autocov_dist(resglm, cbind(dat_heron[,1], dat_heron[,2]),
                  nbs = 50, type = "one", zero.policy = T,
                  style = "B", longlat=TRUE)
dat_heron$RAC<-RAC

m3<-glm(PA~chl_log+smt_100+yft_juv_03_ave+
          bty+bet_juv_03_ave+RAC,
        data=dat_heron, family="binomial")

summary(m3)
print(sum((resid(m3, type="pearson")^2))/df.residual(m3))
print(roc.area(dat_heron$PA, fitted(m3))$A)
pR2(m3)



corm2 <- spline.correlog(dat_heron$Longitude,
        dat_heron$Latitude, residuals(m2,
        type="pearson"), na.rm=T,latlon=T,resamp=10)

corm3 <- spline.correlog(dat_heron$Longitude,
          dat_heron$Latitude, residuals(m3,
          type="pearson"), na.rm=T,latlon=T,resamp=10)


# get odds

exp(cbind(OR = coef(m3), confint(m3)))

#Dominance analysis

dm<-dm2(m3) # works!
#I'll be using McFaddens' r2
dm$contribution.average$r2.m

#make sure it lines up
sum(dm$contribution.average$r2.m)
pR2(m3) # yep :)

#Save dataset we're happy with.
#write.csv(dat_heron, "spreads/dat_heron_3000PA_final.csv", row.names=F, quote=F)

# looking at variable importance/contribution

anova(m3) # anova tests terms sequentially
regTermTest(m3, "chl_log")
anova(m3, update(m3, ~.- "chl_log")) # basically tells the same as anova but gives significance too
Anova(m3) # Anova from car tests terms according to marginality
# i.e. after all other terms have been included

# See if step methods want to shave off or add variables to m3 

for.aic <- step(glm(PA~1,data=dat_heron, family="binomial"),
                direction = "forward", scope = formula(m3), k = 2, trace = 1) # forward AIC
for.bic <- step(glm(PA~1,data=dat_heron, family="binomial"),
                direction = "forward", scope = formula(m3), k = log(nrow(dat_heron)), trace = 0) # forward BIC
back.aic <- step(m3, direction = "backward", k = 2, trace = 0) # backward AIC
back.bic <- step(m3, direction = "backward", k = log(nrow(dat_heron)), trace = 0) # backward BIC
#back steps give glm.fit warning (removing wind?)

formula(for.aic);formula(for.bic);formula(back.aic);formula(back.bic);
# aic methods say don't drop chl, its being weird as its correlated with sst
formula(m2)
# so basically m2 and aic methods are the same
pR2(m3)
pR2(for.aic)
pR2(for.bic) # no diff really in r2

anova(m3, for.aic, for.bic) # same model

# looking at variable importance/contribution
varImp(m3) 
summary(m3)
# seemingly ok, but this is a very crude metric of importance, just scaled
# coefficient and standard error (z value). Collinearity is messing this up.


##### LORD HOWE #####

lhiP<-dat[dat$dset=="LHI" & dat$Count>0,]
lhiA<-dat[dat$dset=="LHI" & dat$Count==0,]
lhiPA<-lhiA[sample(1:nrow(lhiA), 3200),]
dat_lhi<-rbind(lhiP, lhiPA)


# finding distance to nearest point after na.omit, this will identify isolated points
# proj dat into m units for gDistance
sp1<-SpatialPointsDataFrame(dat_lhi[,1:2], data=dat_lhi, proj4string=CRS("+proj=longlat + datum=wgs84"))
sp1proj <- spTransform(sp1, CRS=CRS("+proj=laea +lon_0=155 +lat_0=-21"))

d <- gDistance(sp1proj, byid=T)

min.d <- apply(d, 1, function(x) min(x[x>1000])) # pulls nearest point, which isnt itself or within 1000m projection error range
dat_lhi$min.d<-min.d
hist(dat_lhi$min.d)

plot(Latitude~Longitude, dat_lhi)
points(Latitude~Longitude, dat_lhi[dat_lhi$min.d>50000,], col=2)

dat_lhi<-dat_lhi[dat_lhi$min.d<50000,] # remove orphaned points
# now have a look at collinearity

#pairs(dat_lhi[,-(1:2)], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_lhi[,5:23], upper.panel = panel.smooth,lower.panel=panel.cor)

# final one for model
pairs(dat_lhi[,c(5,6,9,10,12,13,14,16,17, 19,20, 23, 24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# looking at it in full, i would be tempted to remove shd, tmc and yellowfin adult

# Some serious collinearity need to investigate using PCA
#then z transform to standardize
enviro_std<-decostand(dat_lhi[,c(5,6,9,10,12,13,14,16,17, 19,20, 23, 24)], method="standardize")
# takes all varibs but eke also juv and adu tnua forms

# then do pca (just a scaled RDA)
enviro_rda<-rda(enviro_std, scale=T)
summary(enviro_rda, display=NULL)
screeplot(enviro_rda) # badly scaled
#full summary
summary(enviro_rda)

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

# Basically for Heron we have ~50% of the variation explained by PC1 
# Which is effectively a gradient of Latitude with increasing SST and SKJ
# and YFT numbers (Coral Sea) at one end and High Chlorophyll, SHD 
# and deep thermocline at the other end (Tasman sea)

summary(enviro_rda)$species


# check final correltaions and VIFs
pairs(dat_lhi[,c(9,12,13,14,16,17,19 ,23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)

pairs(dat_lhi[,c(9,12,13 ,23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# 
pairs(dat_lhi[,c(5,9,12,19,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
# 

# 
vif(lm(1:nrow(dat_lhi)~bet_adu_03_ave+
         bty+smt_100+wnd+chl_log,
       data=dat_lhi))

# Now we see which variables are better suited to a second deg polynomial

d1<-melt(dat_lhi, id.vars=c("PA", "Count", "dset"))
# Very raw view of how a poly might better suit data than linear
g1<-ggplot(data=d1, aes(y=PA, x=value))
g1+geom_jitter(height=0.1, size=0.5)+geom_smooth(method="glm", colour=2)+
  geom_smooth(method="glm", formula=y~poly(x,2), colour=3)+facet_wrap(~variable, scale="free")

# We would do bl_cw for each variable individually to see the suitability of a poly

m_lin<-glm(PA~bet_adu_03_ave,data=dat_lhi, family="binomial")
m_pol<-glm(PA~poly(bet_adu_03_ave, 2),data=dat_lhi, family="binomial")
d1<-data.frame(PA=dat_lhi$PA, bet_adu_03_ave=dat_lhi$bet_adu_03_ave, m_lin=fitted(m_lin), m_pol=fitted(m_pol))
g1<-ggplot(data=d1, aes(y=PA, x=bet_adu_03_ave))
g1+geom_jitter(height=0.1, size=0.5)+geom_line(aes(y=m_lin, x=bet_adu_03_ave), colour=2)+
  geom_line(aes(y=m_pol, x=bet_adu_03_ave), colour=3)
anova(m_lin, m_pol)
#poly model better in the case of skj_juv it gives the 0 to 1 warning so only use linear
# It was important to work out what was causing the warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
# it was the  skj_juv poly so its now removed and we can proceed!
# not using ekm as its crap anyway

# To account for spatial autocorrelation we calculate autocovariate terms
# using the residuls of the fitted model and a neighbourhood distance to 
# link spatially correlated variables, this term (RAC) is then included in
# the final model to soak up the residual deviance 

# Choose either sst, chla, skj_adu as tropical tuna metric

pairs(dat_lhi[,c(5,9,12,13,16 ,23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)

pairs(dat_lhi[,c(5,9,10,12,13,21 ,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_lhi[,c(5,9,10,12,13 ,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)
pairs(dat_lhi[,c(9,10,12,13, 23,24) ], upper.panel = panel.smooth,lower.panel=panel.cor)


# check which tuna metric best predictor

AIC(glm(PA~skj_adu_03_ave,data=dat_lhi, family="binomial"),
    glm(PA~skj_juv_03_ave,data=dat_lhi, family="binomial"),
    glm(PA~bet_adu_03_ave,data=dat_lhi, family="binomial"),
    glm(PA~bet_juv_03_ave,data=dat_lhi, family="binomial"),
    glm(PA~yft_adu_03_ave,data=dat_lhi, family="binomial"),
    glm(PA~yft_juv_03_ave,data=dat_lhi, family="binomial"))

# bet_adu, bet_juv, skj_juv

# see if there are outliers in selected variables
qplot(data=dat_lhi, x=skj_juv_03_ave, bins=50)+facet_grid(PA~.)

# chl_log not important very, get pushed around by tuna metrics

m2a<-glm(PA~chl_log+smt_100+bty+bet_juv_03_ave,
         data=dat_lhi, family="binomial")

#  bty and yft are correlated but seem ok together
m2b<-glm(PA~smt_100+bty+bet_juv_03_ave+yft_adu_03_ave,
         data=dat_lhi, family="binomial")

m2c<-glm(PA~smt_100+bty+bet_adu_03_ave,
         data=dat_lhi, family="binomial")

anova(m2a, m2b, m2c);AIC(m2a, m2b, m2c)
pR2(m2a)[4];pR2(m2b)[4];pR2(m2c)[4]

m2<-glm(PA~smt_100+bty+bet_juv_03_ave+yft_adu_03_ave,
         data=dat_lhi, family="binomial");resglm<-residuals(m2, type="pearson")

summary(m2)

print(sum((resid(m2, type="pearson")^2))/df.residual(m2))
print(roc.area(dat_lhi$PA, fitted(m2))$A)
pR2(m2)

resglm<-residuals(m2, type="pearson")
sp2<-SpatialPointsDataFrame(dat_lhi[,1:2], data=data.frame(resglm), 
                            proj4string=CRS("+proj=longlat + datum=wgs84"))
bubble(sp2, zcol='resglm')
#bubble(sp2, zcol='resglm', do.sqrt=F)

#write.csv(cbind(dat_lhi, resglm), "test5.csv", quote=FALSE)


#corglm <- correlog(dat_lhi$Longitude, dat_lhi$Latitude, residuals(m2, type="pearson"), na.rm=T,
#                   latlon=T, increment=25,resamp=1)
RAC<-autocov_dist(resglm, cbind(dat_lhi[,1], dat_lhi[,2]),
                  nbs = 50, type = "one", zero.policy = T,
                  style = "B", longlat=TRUE)

dat_lhi$RAC<-RAC

m3<-glm(PA~smt_100+bty+
          bet_juv_03_ave+yft_adu_03_ave+RAC,
        data=dat_lhi, family="binomial")

plot(m3)
print(sum((resid(m3, type="pearson")^2))/df.residual(m3))

# not too bad overdisp but still some big negative resids

dat_lhi<-dat_lhi[-which(resid(m3, type="pearson")< -10),]
#refit the model

m2<-glm(PA~smt_100+bty+bet_juv_03_ave+yft_adu_03_ave,
        data=dat_lhi, family="binomial");resglm<-residuals(m2, type="pearson")

RAC<-autocov_dist(resglm, cbind(dat_lhi[,1], dat_lhi[,2]),
                  nbs = 50, type = "one", zero.policy = T,
                  style = "B", longlat=TRUE)

dat_lhi$RAC<-RAC

m3<-glm(PA~smt_100+bty+
          bet_juv_03_ave+yft_adu_03_ave+RAC,
        data=dat_lhi, family="binomial")

summary(m3)
print(sum((resid(m3, type="pearson")^2))/df.residual(m3))
print(roc.area(dat_lhi$PA, fitted(m3))$A)
pR2(m3)


corm2 <- spline.correlog(dat_lhi$Longitude, dat_lhi$Latitude, residuals(m2, type="pearson"), na.rm=T,
                  latlon=T,resamp=10)
corm3 <- spline.correlog(dat_lhi$Longitude, dat_lhi$Latitude, residuals(m3, type="pearson"), na.rm=T,
                  latlon=T,resamp=10)

# get odds

exp(cbind(OR = coef(m3), confint(m3)))

# run dominanace analysis

m3d<-dm2(m3)
m3d$contribution.average$r2.m

# Save out spac plots

out_spac<-rbind(
data.frame(dset="heron_m2", Dist=corm2$boot$boot.summary$predicted$x[1,],
           SPAC_025=corm2$boot$boot.summary$predicted$y[3,],
           SPAC_Ave=corm2$boot$boot.summary$predicted$y[6,],
           SPAC_95=corm2$boot$boot.summary$predicted$y[9,]),
data.frame(dset="heron_m3", Dist=corm3$boot$boot.summary$predicted$x[1,],
           SPAC_025=corm3$boot$boot.summary$predicted$y[3,],
           SPAC_Ave=corm3$boot$boot.summary$predicted$y[6,],
           SPAC_95=corm3$boot$boot.summary$predicted$y[9,]),
data.frame(dset="lhi_m2", Dist=corm4$boot$boot.summary$predicted$x[1,],
           SPAC_025=corm4$boot$boot.summary$predicted$y[3,],
           SPAC_Ave=corm4$boot$boot.summary$predicted$y[6,],
           SPAC_95=corm4$boot$boot.summary$predicted$y[9,]),
data.frame(dset="lhi_m3", Dist=corm5$boot$boot.summary$predicted$x[1,],
           SPAC_025=corm5$boot$boot.summary$predicted$y[3,],
           SPAC_Ave=corm5$boot$boot.summary$predicted$y[6,],
           SPAC_95=corm5$boot$boot.summary$predicted$y[9,]))

qplot(data=out_spac, x=Dist, y=SPAC_Ave, colour=dset, geom="line")+
      geom_ribbon(aes(ymin=SPAC_025, ymax=SPAC_95, fill=dset),alpha=0.25)+theme_classic()
#sweet
#write.csv(out_spac, "paper_results/climatology_SPAC_results.csv", quote=F, row.names=F)

##### Making plots #####

# Make boxplots comparing oceanic and tuna covariates between
# years and colonies.

datv4<-read.csv("spreads/paper2_extractionV4.csv", h=T, strip.white=T)
# use v4 as has individual years split

pres_all<-datv4[datv4$dtyp=="ud50_pres",]




#! Understanding outputs !#

# make same final model here but with non orthanonol polynomial,
# this means that the coefficients are scaled correctly and more interpretable.

m3_interp<-glm(PA~chl_log+wnd+
                 smt_100+poly(bty,2, raw=T)+
                 bet_adu_03_ave+
                 +sst+RAC,
               data=dat_heron, family="binomial")
anova(m3, m3_interp)

summary(m3_interp)



-53.26117519/(-28.50752856*2)
# 0.9341598
-1.922579e-03/(-2.427522e-07*2)
# 3959.962

-226.254389/(-86.476828*2)

# Working on my understanding of interpreting odds ratio

wnd<-c(4,4.2,4.3,4.2,4.8,4.3,5,5.3,5.4,5.6,5.8,5.9,6, 6.1,6.4)
PA<-c(0,0,0,0,1,0,0,1,1,0,1,1,1,1,1)

testdat<-data.frame(PA=PA,wnd=wnd)
summary(glm.mod <- glm(PA~ wnd, family = "binomial", data=testdat))

## logit = -0.8690 + (-1.0769) * x
v1 <- -16.611  + (3.268)*5 # log odds of Presence at wind=5
v2 <- -16.611  + (3.268)*6 # log odds of Presence at wind=6

v1
# -0.271
v2
# 2.997
v2-v1
# 3.268 ok so 1 unit change in wind increases log odds of
# Presence by 3.268 - exactly what model summary says!
# 

# odds of v1 and v2 are
exp(v1)
# 0.7626165
exp(v2)
# 20.02537

exp(v1)/(1+exp(v1))
# 0.4326616 Probability of presence at wnd=5
exp(v2)/(1+exp(v2))
# 0.9524384 Probability of presence at wnd=6
         
exp(coef(m3))
plogis(coef(m3))

# Loop to check model predictions are correct, this checks for collinearity
# to see if some trends look wrong or if other variables are seeing their response
# supressed due to collinearity with others. For presentation we could show model m2
# Where the RAC term is not included and therefore the responses of meaningful variables
# are clearer.

names(model.frame(m3)) # running on m2!
for (i in c("chl_log", "sst","wnd","tmc","smt_100",               
            "bty","bet_adu_03_ave","RAC"))
{
  temp<-dat_lhi[1,-4]
  temp[1,]<-as.vector(apply(dat_lhi[,-4], 2,median))
  
  heron_pred<-data.frame(varib=dat_lhi[,which(names(dat_lhi)==i)], 
                         temp)
  names(heron_pred)[names(heron_pred)==i]<-"nope"
  names(heron_pred)[names(heron_pred)=="varib"]<-i
  
  p1<-predict(m2, newdata=heron_pred, type="response")
  
  d1<-data.frame(PA=dat_lhi$PA, varib=heron_pred[,1], pred=p1)
  g1<-ggplot(data=d1, aes(y=PA, x=varib))
  print(i)
  print(g1+geom_jitter(height=0.1, size=0.5)+
          geom_line(aes(y=pred, x=varib), colour=2))
  readline("hi@")
}




####### ~~~ Unused ~~~~ #######

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

me.fit<-ME(PA~shd+ekm+chl_log+wnd+tmc+smt_100+bty+
             bet_juv_03_ave+yft_adu_03_ave+
             skj_juv_03_ave, data=dat_heron, 
           family=binomial, listw = dat.wt4, verbose=T,alpha=0.05 )


