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

env_vars<-c(sst, shd, ekm, chl, wnd, tmc, smt, bty)

#read in seapodym data

bet_adu<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/bet_adu_03_ave.tif")
bet_juv<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/bet_juv_03_ave.tif")
bet_tot<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/bet_tot_03_ave.tif")
skj_adu<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/skj_adu_03_ave.tif")
skj_juv<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/skj_juv_03_ave.tif")
skj_tot<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/skj_tot_03_ave.tif")
yft_adu<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/yft_adu_03_ave.tif")
yft_juv<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/yft_juv_03_ave.tif")
yft_tot<-raster("/home/mark/grive/phd/sourced_data/SEAPODYM/month_ave_intermin_1deg_ref2015/yft_tot_03_ave.tif")

tun_stack<-stack(bet_adu, bet_juv, bet_tot, skj_adu, skj_juv, skj_tot, yft_adu, yft_juv, yft_tot)

kernels<-list.files("spatial")
kernels<-kernels[grep("shp", kernels)]
kernels<-kernels[-c(9,10)]

#Construct standardised sampling grid so even and fair of sampling of different resolution datasets
# set resolution at 0.25 degree, good comprimise between 0.01 bathy and 1 deg tuna or, 2 deg thermo.. hmm not sure about that one

ext_all<-NULL
for( i in kernels)
  {
  sp1<-readOGR(dsn=paste("/home/mark/grive/phd/analyses/paper2/spatial/",
              i, sep=""), layer=substr(i, 1, (nchar(i)-4)))
  
  if(length(grep("LT", i))==1) # grep throws a integer(0), use length() to capture in if statement
    {
    r1<-rasterize(sp1[3,], chl) #makes raster of 50% UD at 0.1 deg resolution
    p1<-rasterToPoints(r1)
    dtp="ud50_pres"
    }else{
    r1<-rasterize(spsample(sp1, n=6000, type="random"), chl) # using rougly 3:1 background sample.. still zero inf?
    p1<-rasterToPoints(r1)
    dtp="psuedo_abs"
    }
        
    #plot(sp1[3,])
    #plot(SpatialPoints(p1), add=T, col=2)
  ext_v<-data.frame(Longitude=p1[,1], Latitude=p1[,2], dtyp=dtp, dset=substr(i, 1, (nchar(i)-4)), 
                      sst=1, shd=1, ekm=1, chl=1, wnd=1, tmc=1, smt=1, bty=1)
  counter=0
  for(j in env_vars )
      {
      counter<-counter+1
      e1<-extract(j, ext_v[,1:2]) # extract all pixels within 20 km of point and average
      ext_v[,4+counter]<-e1
      }
     
  e2<-extract(tun_stack, ext_v[,1:2]) # no point using the buffer here as pixels at 1 degree
  ext_v<-cbind(ext_v, e2)
  ext_all<-rbind(ext_all, ext_v)  
}   

write.csv(ext_all, "spreads/paper2_extractionV2.csv", quote=F, row.names=F)
#writing out data, note as ive randomly sampled the psuedo abs data, it might be
# a good idea to iterate modelling with additional random samples of more or less
# to see when results stabalise... ie we have sampled the area comprehensivly







