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

dat<-read.csv("~/grive/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm_attribs.csv", h=T)

# need to define h value (Scale), thinking max dist in data between 2 points (10 mins)
max(dat$Vel)/6 #max velocity is 74 km/h, so is 12km in 10 mins (resmaple time)
hscale<-12

source("~/grive/phd/scripts/github_MIBA/batchUD.R")

#source("~/grive/phd/scripts/MIBA_scripots_revised/BatchUD_revised.r")

for(j in c(99,75,50,25))
{
  dat$ID<-j
  UD_out<-batchUD(dat[with(dat, trip_type=="L" &
                  Colony=='Heron' & Year=="2015"),],
                  Scale = hscale, UDLev = j)
  
  if(j==99){all_UDs<-UD_out}else{all_UDs<-spRbind(all_UDs, UD_out)}
  
  plot(all_UDs, border=factor(all_UDs$id))
}

all_UDs <- spTransform(all_UDs, CRS=CRS("+proj=longlat +ellps=WGS84"))
#I updated Lubuntu and now proj is handed correctly so can deal with WGS!! :)

plot(all_UDs, border=factor(all_UDs$id), lwd=2)

writeOGR(all_UDs, layer="trial", dsn="spatial", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)

## to define area available to each colony constuct a radius (using max trip dist)
## then clip out land

d1<-max(dat$ColDist) #1183621 m

Colony_heron<-SpatialPoints(data.frame(Longitude=151.913773, Latitude=-23.44306 ), proj4string=CRS("+proj=longlat + datum=wgs84"))
Colony_lhi<-SpatialPoints(data.frame(Longitude=159.05991, Latitude=-31.52459 ), proj4string=CRS("+proj=longlat + datum=wgs84"))

DgProj <- CRS("+proj=laea +lon_0=156 +lat_0=-17")

heronProj <- spTransform(Colony_heron, CRS=DgProj)
lhiProj <- spTransform(Colony_lhi, CRS=DgProj)

#need to use a width as specified above
heronBuffProj <- gBuffer(heronProj, width=d1, quadsegs=50)
lhiBuffProj <- gBuffer(lhiProj, width=d1, quadsegs=50)
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
plot(heronBuffWgsclip)

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

heron_psuedo<-spsample(heronBuffWgsclip, 2000, type="random")
plot(heron_psuedo)
lhi_psuedo<-spsample(lhiBuffWgsclip, 2000, type="random")

#read in climatology env data

sst<-raster("~/grive/phd/sourced_data/env_data/climatologies/CPCsstn2016-03-01.nc")
shd<-raster("~/grive/phd/sourced_data/env_data/climatologies/CTCsshd2016-03-01.nc")
ekm<-raster("~/grive/phd/sourced_data/env_data/climatologies/CQCwekm2016-03-01.nc")
chl<-raster("~/grive/phd/sourced_data/env_data/climatologies/CPMbfp12016-03-01.nc")
wnd<-raster("~/grive/phd/sourced_data/env_data/climatologies/TNCwwww1995-03-31.nc")
tmc<-raster("~/grive/phd/sourced_data/env_data/climatologies/thermocline_3.tif")
smt<-raster("/home/mark/grive/phd/sourced_data/env_data/seamounts/d_seamounts_wgs.tif")
bty<-raster("~/grive/phd/sourced_data/env_data/phd_bathy/GRIDONE_2D_100.0_-45.0_180.0_40.0.nc")




