# 13/09/16 Lancaster, UK
# constuct kernels from GPS data and psuedo-absence foraging range then 
# extract from climatology data

rm(list=ls())

library(sp)
library(maptools)
library(rgdal)
library(adehabitat)                        
library(geosphere) 
library(rgeos)

setwd("~/grive/phd/analyses/paper2")

######### ************** PSUEDO ABSENCE CONSTRUCTION ************** ##########

## to define area available to each colony constuct a radius (using max trip dist)
## then clip out land

Colony_heron<-SpatialPoints(data.frame(Longitude=151.913773, Latitude=-23.44306 ), proj4string=CRS("+proj=longlat + datum=wgs84"))
Colony_lhi<-SpatialPoints(data.frame(Longitude=159.05991, Latitude=-31.52459 ), proj4string=CRS("+proj=longlat + datum=wgs84"))

DgProj <- CRS("+proj=laea +lon_0=156 +lat_0=-17")

heronProj <- spTransform(Colony_heron, CRS=DgProj)
lhiProj <- spTransform(Colony_lhi, CRS=DgProj)

#need to use a width as specified above # now 1400m as Heron 2013 PTT data shows
heronBuffProj <- gBuffer(heronProj, width=1400000, quadsegs=50)
lhiBuffProj <- gBuffer(lhiProj, width=1400000, quadsegs=50)
#TBuffProj@polygons[[1]]@ID <- as.character(i)
  
heronBuffWgs <- spTransform(heronBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))
lhiBuffWgs <- spTransform(lhiBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))
 # now clip out land help here: http://gis.stackexchange.com/questions/109639/how-to-reverse-clip-erase-in-r

# for some reason readOGR doesnt like the "~/grive" dir 
world<-readOGR(dsn="/home/mark/grive/phd/sourced_data/env_data/global_map/ne_110m_land.shp", layer="ne_110m_land")

plot(world)
plot(lhiBuffWgs, add=T)
plot(heronBuffWgs, add=T, border=2)

heronBuffWgsclip<-gDifference(heronBuffWgs, world) # warning but should be ok (bassically the smae proj)
lhiBuffWgsclip<-gDifference(lhiBuffWgs, world) # warning but should be ok (bassically the smae proj)
plot(heronBuffWgsclip) # edited out bit beyond cape york in qgis

## extra clip out bottom of lhi range against birdlife range map
notrange<-readOGR(dsn="/home/mark/grive/phd/sourced_data/otherGIS/not_WTSH_range.shp", layer="not_WTSH_range")

lhiBuffWgsclip<-gDifference(lhiBuffWgsclip, notrange) 
plot(lhiBuffWgsclip)

writeOGR(SpatialPolygonsDataFrame(heronBuffWgsclip, data.frame(ID=1)), layer="heronBuff", dsn="spatial", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)
writeOGR(SpatialPolygonsDataFrame(lhiBuffWgsclip, data.frame(ID=1)), layer="lhiBuff", dsn="spatial", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)

# coolio now set about extraction, let raster do random point sample over access from colony weight
library(raster)
# what we're gonna do is randomly sample over an access gradient. Figured this idea while
# at birdlife basically we use basic sample and enter a distance cost raster around each colony
# as a dataframe and let sample do the rest!

# we need as fine-a resolution template raster
bathy<-raster("~/grive/phd/sourced_data/env_data/phd_bathy/GRIDONE_2D_100.0_-45.0_180.0_40.0.nc")

heronRas<-crop(bathy, extent(heronBuffWgsclip))
heronRas<-rasterize(Colony_heron, heronRas)
heronD<-distance(heronRas)
heronD[heronD==0,]<-1702 # removes 0 distance pixel and sets it to min pixel distance ~ 2000m
heronD<-1/heronD # reverses distance surface to make access surface
heronApts<-rasterToPoints(heronD, spatial=T)
heronApts@ proj4string<-CRS(proj4string(heronBuffWgsclip)) #force to same projection (they are already)
heronApts<-heronApts[heronBuffWgsclip,] # this is a clip!! thanks to for such simple solution http://robinlovelace.net/r/2014/07/29/clipping-with-r.html

plot(heronApts[sample(length(heronApts),2000, prob=heronApts$layer),]) # layer is our access surface

### after all this I've actually thought myself out of using the weighting!!
### going to just randomly distribute points over the entire area.

#heron_psuedo<-spsample(heronBuffWgsclip, 2000, type="random")
#plot(heron_psuedo)
#lhi_psuedo<-spsample(lhiBuffWgsclip, 2000, type="random")

######### ************** KERNEL CONSTRUCTION ************** ##########

dat<-read.csv("~/grive/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm_attribs.csv", h=T)
heron_dat<-read.csv("~/grive/phd/sourced_data/Heron/heron_paper_Scarla/fi_data_04_sep/All LT PTT 2012 & 2013 day foraging data for kernels & 2006 2011.csv", h=T)

# this bif of code cuts out a 300km buffer from the PTT data as per mcduie et al
heron_old<-SpatialPointsDataFrame(SpatialPoints(data.frame(heron_dat$Longitude, heron_dat$Latitude),
                                 proj4string=CRS( "+proj=longlat +ellps=WGS84")), data=heron_dat)

heronBuffProj <- gBuffer(heronProj, width=300000, quadsegs=50) # 300km buffer as per mcduie et al
heronBuffWgs <- spTransform(heronBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))
plot(heron_old);plot(heronBuffWgs, add=T, border=2)
ov1<-over(heron_old, heronBuffWgs) # Using over solves lots of problems for spatial subsetting http://gis.stackexchange.com/questions/63793/how-to-overlay-a-polygon-over-spatialpointsdataframe-and-preserving-the-spdf-dat
heron_old<-heron_old[which(is.na((ov1))),]
plot(heron_old)


dat_old<-SpatialPointsDataFrame(SpatialPoints(data.frame(dat$Longitude, dat$Latitude),
                                proj4string=CRS( "+proj=longlat +ellps=WGS84")), data=dat)

datBuffProj <- gBuffer(lhiProj, width=20000, quadsegs=50) # first cut out 20 lhi buffer
datBuffWgs <- spTransform(datBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))
plot(dat_old);plot(datBuffWgs, add=T, border=2)
ov1<-over(dat_old, datBuffWgs) # Using over solves lots of problems for spatial subsetting http://gis.stackexchange.com/questions/63793/how-to-overlay-a-polygon-over-spatialpointsdataframe-and-preserving-the-spdf-dat
dat_old<-dat_old[which(is.na((ov1))),]
plot(dat_old)

datBuffProj <- gBuffer(heronProj, width=20000, quadsegs=50) # now cut out 20km heron buffer
datBuffWgs <- spTransform(datBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))
plot(dat_old);plot(datBuffWgs, add=T, border=2)
ov1<-over(dat_old, datBuffWgs) # Using over solves lots of problems for spatial subsetting http://gis.stackexchange.com/questions/63793/how-to-overlay-a-polygon-over-spatialpointsdataframe-and-preserving-the-spdf-dat
dat_old<-dat_old[which(is.na((ov1))),]
plot(dat_old) # FINALLY!!

## Deciding on H value for kernels
d1<-max(dat$ColDist) #1183621 m
# need to define h value (Scale), thinking max dist in data between 2 points (10 mins)
max(dat$Vel)/6 #max velocity is 74 km/h, so is 12km in 10 mins (resmaple time)
hscale<-20 # following McDuie et al 2015

source("~/grive/phd/scripts/github_MIBA/batchUD.R")

# decided to use 2011, 2013 PTT data, 2006 and 2012 seems to small 
#heron_old[heron_old$Year==2011,]
#dat_old[with(dat_old, trip_type=="L" & Colony=='Heron' & Year=="2015"),]
heron_old<-heron_old@data
# while we have this data we'll work out how long each LT is active for, for trip summary table 2
aggregate(Date~Nest_id, heron_old[heron_old$Year==2011,], FUN =function(x){length(unique(x))})
mean(c(11,11,7)); sd(c(11,11,7))

aggregate(Date~Nest_id, heron_old[heron_old$Year==2013,], FUN =function(x){length(unique(x))})
mean(c(7,9,9,7,15,11,7,4,6)); sd(c(7,9,9,7,15,11,7,4,6))

dat_old<-dat_old@data

for(j in c(99,75,50,25))
{
  #heron_old$ID<-j
  dat_old$ID<-j
  UD_out<-batchUD(dat_old[with(dat_old, trip_type=="L" & Colony=='LHI' & Year=="2016"),],
                  Scale = hscale, UDLev = j)
  
  if(j==99){all_UDs<-UD_out}else{all_UDs<-spRbind(all_UDs, UD_out)}
  
  plot(all_UDs, border=factor(all_UDs$id))
}

all_UDs <- spTransform(all_UDs, CRS=CRS("+proj=longlat +ellps=WGS84"))
#I updated Lubuntu and now proj is handed correctly so can deal with WGS!! :)

plot(all_UDs, border=factor(all_UDs$id), lwd=2)

writeOGR(all_UDs, layer="LTLHIGPS2016", dsn="spatial", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)


#read in climatology env data

sst<-raster("~/grive/phd/sourced_data/env_data/climatologies/CPCsstn2016-03-01.nc")
shd<-raster("~/grive/phd/sourced_data/env_data/climatologies/CTCsshd2016-03-01.nc")
ekm<-raster("~/grive/phd/sourced_data/env_data/climatologies/CQCwekm2016-03-01.nc")
chl<-raster("~/grive/phd/sourced_data/env_data/climatologies/CPMbfp12016-03-01.nc")
wnd<-raster("~/grive/phd/sourced_data/env_data/climatologies/TNCwwww1995-03-31.nc")
tmc<-raster("~/grive/phd/sourced_data/env_data/climatologies/thermocline_3.tif")
smt<-raster("/home/mark/grive/phd/sourced_data/env_data/seamounts/d_seamounts_wgs.tif")
bty<-raster("~/grive/phd/sourced_data/env_data/phd_bathy/GRIDONE_2D_100.0_-45.0_180.0_40.0.nc")

tmc<-resample(tmc, shd, method='bilinear') # resmaple tmc to shd resolution .25

env_vars<-c(sst, shd, ekm, chl, wnd, tmc, smt, bty)

#read in seapodym data

bet_adu<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/bet_adu_03_ave.tif"),
                  shd, method='bilinear') # resmaple tmc to shd resolution .25
bet_juv<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/bet_juv_03_ave.tif"),
                  shd, method='bilinear')
bet_tot<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/bet_tot_03_ave.tif"),
                  shd, method='bilinear')
skj_adu<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/skj_adu_03_ave.tif"),
                  shd, method='bilinear')
skj_juv<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/skj_juv_03_ave.tif"),
                  shd, method='bilinear')
skj_tot<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/skj_tot_03_ave.tif"),
                  shd, method='bilinear')
yft_adu<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/yft_adu_03_ave.tif"),
                  shd, method='bilinear')
yft_juv<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/yft_juv_03_ave.tif"),
                  shd, method='bilinear')
yft_tot<-resample(raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/yft_tot_03_ave.tif"),
                  shd, method='bilinear')

tun_stack<-stack(bet_adu, bet_juv, bet_tot, skj_adu, skj_juv, skj_tot, yft_adu, yft_juv, yft_tot)

kernels<-list.files("spatial")
kernels<-kernels[grep("shp", kernels)]
kernels<-kernels[-c(9,10)]

#Construct standardised sampling grid so even and fair of sampling of different resolution datasets
# set resolution at 0.25 degree, good comprimise between 0.01 bathy and 1 deg tuna or, 2 deg thermo.. hmm not sure about that one

for( i in kernels[grep("LTLHI", kernels)])
  {
  sp1<-readOGR(dsn=paste("/home/mark/grive/phd/analyses/paper2/spatial/",
              i, sep=""), layer=substr(i, 1, (nchar(i)-4)))
  
  r1<-rasterize(sp1[3,], chl, field=1, background=0) 
  if(which(kernels[grep("LTLHI", kernels)]==i)==1){
    rsum=r1}else{rsum=rsum+r1}
  }
p1<-rasterToPoints(rsum, fun=function(x){x>0})

lhiP<-data.frame(Longitude=p1[,1], Latitude=p1[,2], Count=p1[,3], dset="LHI", 
                  sst=1, shd=1, ekm=1, chl=1, wnd=1, tmc=1, smt=1, bty=1)

for( i in kernels[grep("LTHeron", kernels)])
{
  sp1<-readOGR(dsn=paste("/home/mark/grive/phd/analyses/paper2/spatial/",
                         i, sep=""), layer=substr(i, 1, (nchar(i)-4)))
  
  r1<-rasterize(sp1[3,], chl, field=1, background=0) 
  if(which(kernels[grep("LTHeron", kernels)]==i)==1){
    rsum=r1}else{rsum=rsum+r1}
}
p1<-rasterToPoints(rsum, fun=function(x){x>0})

herP<-data.frame(Longitude=p1[,1], Latitude=p1[,2], Count=p1[,3], dset="Heron", 
                 sst=1, shd=1, ekm=1, chl=1, wnd=1, tmc=1, smt=1, bty=1)

sp1<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/lhiBuff.shp",
                   layer="lhiBuff")

r1<-rasterize(sp1, chl) # here we extract all points and will do sample() testing as part of the modelling
p1<-rasterToPoints(r1)

lhiA<-data.frame(Longitude=p1[,1], Latitude=p1[,2], Count=0, dset="LHI", 
                 sst=1, shd=1, ekm=1, chl=1, wnd=1, tmc=1, smt=1, bty=1)

sp1<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/heronBuff.shp",
             layer="heronBuff")

r1<-rasterize(sp1, chl) # here we extract all points and will do sample() testing as part of the modelling
p1<-rasterToPoints(r1)

herA<-data.frame(Longitude=p1[,1], Latitude=p1[,2], Count=0, dset="Heron", 
                 sst=1, shd=1, ekm=1, chl=1, wnd=1, tmc=1, smt=1, bty=1)

 ext_v<-rbind(lhiP, lhiA, herP, herA)       
    #plot(sp1[3,])
    #plot(SpatialPoints(p1), add=T, col=2)
  counter=0
  for(j in env_vars )
      {
      counter<-counter+1
      e1<-extract(j, ext_v[,1:2]) # extract all pixels within 20 km of point and average
      ext_v[,4+counter]<-e1
      }
     
  e2<-extract(tun_stack, ext_v[,1:2]) # no point using the buffer here as pixels at 1 degree
  ext_v<-cbind(ext_v, e2)
 

write.csv(ext_v, "spreads/paper2_extractionV5.csv", quote=F, row.names=F)
#writing out data, note as ive randomly sampled the psuedo abs data, it might be
# a good idea to iterate modelling with additional random samples of more or less
# to see when results stabalise... ie we have sampled the area comprehensivly

# Make spatial lines of all tracking data for paper figure map

dat<-read.csv("~/grive/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm_attribs.csv", h=T)
heron_dat<-read.csv("~/grive/phd/sourced_data/Heron/heron_paper_Scarla/fi_data_04_sep/All LT PTT 2012 & 2013 day foraging data for kernels & 2006 2011.csv", h=T)

d1<-data.frame(Latitude=dat[dat$trip_type=="L",]$Latitude, Longitude=dat[dat$trip_type=="L",]$Longitude,
               trip_id=as.character(dat[dat$trip_type=="L",]$trip_id), Year=dat[dat$trip_type=="L",]$Year,
               Colony=dat[dat$trip_type=="L",]$Colony)

d2<-data.frame(Latitude=heron_dat$Latitude, Longitude=heron_dat$Longitude,
               trip_id=as.character(heron_dat$Nest_id), Year=heron_dat$Year,
               Colony="Heron")

d2<-d2[d2$Year=="2011" | d2$Year=="2013",]

d3<-rbind(d1, d2)

#get centroids of each dataset
d3$ID<-paste(d3$Colony, d3$Year)
aggregate(Latitude~ID, data=d3, FUN=mean)
aggregate(Longitude~ID, data=d3, FUN=mean)

ID  Latitude
#1 Heron 2011 -21.1, 154.7
#2 Heron 2013 -20.9, 153.4
#3 Heron 2015 -21.0, 156.0
#4   LHI 2014 -30.9, 156.5
#5   LHI 2015 -27.7, 158.2
#6   LHI 2016 -31.9, 157.3


trip_IDs<-unique(d3$trip_id)

year_id<-NULL
colony_id<-NULL
for(i in trip_IDs)
{
  Trip <- d3[d3$trip_id == i,]
  L1 <- Line(as.matrix(data.frame(Trip$Longitude,Trip$Latitude)))
  Ls1 <- Lines(L1, ID=which(trip_IDs == i))
  SpLs1 <- SpatialLines(list(Ls1), CRS("+proj=longlat"))
  if(which(trip_IDs == i) == 1) {SpLZ <- SpLs1} else
    SpLZ <- spRbind(SpLZ, SpLs1)
  year_id<-c(year_id, unique(Trip$Year))
  colony_id<-c(colony_id, unique(Trip$Colony))
}


Tbl <- data.frame(Name_0 = 1, Name_1 = 1:length(trip_IDs), ID = trip_IDs, Year = year_id, Colony= colony_id)
row.names(Tbl) <- Tbl$Name_1
SLDF<-SpatialLinesDataFrame(SpLZ, Tbl)  

setwd("~/grive/phd/analyses/paper2")

writeOGR(SLDF, layer="paper_tracking_lines", dsn="spatial", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)


# areas of 99% UD: earch area

k1<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2014.shp",
            layer="LTLHIGPS2014")

k2<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2015.shp",
            layer="LTLHIGPS2015")

k3<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2016.shp",
            layer="LTLHIGPS2016")

h1<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTHeronPTT2011.shp",
            layer="LTHeronPTT2011")

h2<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTHeronPTT2013.shp",
            layer="LTHeronPTT2013")

h3<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTHeronGPS2015.shp",
            layer="LTHeronGPS2015")


# get centroids of 50% UDs

rbind(data.frame(ID="LHI14", gCentroid(k1[3,])),
      data.frame(ID="LHI15", gCentroid(k2[3,])),
      data.frame(ID="LHI16", gCentroid(k3[3,])),
      data.frame(ID="HER11", gCentroid(h1[3,])),
      data.frame(ID="HER13", gCentroid(h2[3,])),
      data.frame(ID="HER15", gCentroid(h3[3,])))

#ID        x         y
#1  LHI14 156.1532 -31.83339
#11 LHI15 158.0625 -27.72413
#12 LHI16 157.1900 -31.83591
#13 HER11 155.7414 -19.20367
#14 HER13 154.2389 -19.78477
#15 HER15 156.9149 -20.40979



library(raster)
# function area can calc to m2 from lat long

area(k1[1,])/1000000
area(k2[1,])/1000000
area(k3[1,])/1000000
area(h1[1,])/1000000
area(h2[1,])/1000000
area(h3[1,])/1000000





## attempt at repeatability of kernels

k1<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2014.shp",
            layer="LTLHIGPS2014")

k2<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2015.shp",
            layer="LTLHIGPS2015")

k3<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2016.shp",
            layer="LTLHIGPS2016")

k1<-k1[3,]
k1@polygons[[1]]@ID<-"2015" ; row.names(k1@data)<-"2015"
k2<-k2[3,]
k2@polygons[[1]]@ID<-"2013" ; row.names(k2@data)<-"2013"
k3<-k3[3,]
k3@polygons[[1]]@ID<-"2011" ; row.names(k3@data)<-"2011"



all_UDs<-spRbind(k1, k2)
all_UDs<-spRbind(all_UDs, k3)
all_UDs$id<-c("2015", "2013", "2011")
plot(all_UDs, border=factor(all_UDs$id))


heron_pols<-disaggregate(all_UDs)

DgProj <- CRS("+proj=laea +lon_0=155 +lat_0=-22")

heron_pols <- spTransform(heron_pols, CRS=DgProj)

library(reshape2)
haw_test<-rbind(
data.frame(year="2011-2013",Dists=melt(gDistance(heron_pols[heron_pols$id==2011,],
                 heron_pols[heron_pols$id==2013,], 
                 hausdorff=TRUE, byid=TRUE))[,3]),
data.frame(year="2013-2015",Dists=melt(gDistance(heron_pols[heron_pols$id==2013,],
                 heron_pols[heron_pols$id==2015,], 
                 hausdorff=TRUE, byid=TRUE))[,3]),
data.frame(year="2011-2015",Dists=melt(gDistance(heron_pols[heron_pols$id==2011,],
                 heron_pols[heron_pols$id==2015,], 
                 hausdorff=TRUE, byid=TRUE))[,3]))

hist(haw_test$Dists)
#hist(sqrt(haw_test$Dists))
qplot(data=haw_test, x=year, y=Dists, geom="boxplot")

haw_test$year<-relevel(haw_test$year, ref="2011-2015")
m1<-lm(Dists~year, haw_test)
summary(m1)


source("~/grive/phd/scripts/MIBA_scripts_revised/varianceTest_revised.r")
bird_string<-as.character(heron_pols$id)

bird_string<-as.character(seq(1:length(heron_pols)))

vt<-varianceTest(heron_pols, bird_string, Iteration=10)
vt
#hmm need to figue out this one

heron_pols


#Or some form of repeatabililty or general statisitcal test to see difference between years



