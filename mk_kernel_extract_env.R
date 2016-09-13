# 13/09/16 Lancaster, UK
# constuct kernels from GPS data and psuedo-absence foraging range then 
# extract from climatology data

source("~/grive/phd/scripts/MIBA_scripts_revised/BatchUD.r")

# need to define h value (Scale), thinking max dist in data between 2 points (10 mins)

for(j in c(99,75,50,25))
{
  final_out$ID<-j
  UD_out<-batchUD(na.omit(final_out[with(final_out, trip_type=="S" &
                                           Colony=='Heron' & Year=="2015"),]),
                  Scale = S_scale, UDLev = j)
  
  if(j==99){all_UDs<-UD_out}else{all_UDs<-spRbind(all_UDs, UD_out)}
  
  plot(all_UDs, border=factor(all_UDs$ID))
}

all_UDs <- spTransform(all_UDs, CRS=CRS("+proj=longlat +ellps=WGS84"))

plot(all_UDs, border=factor(all_UDs$ID), lwd=2)
map("world", add=T, fill=T, col="darkolivegreen3")

writeOGR(all_UDs, layer="Heron_2015_ST", dsn="D:/research/phd/results", driver="ESRI Shapefile", verbose=TRUE, overwrite=T)

## to define area available to each colony constuct a radius (using max trip dist)
## then clip out land

Colony_heron<-SpatialPoints(data.frame(Longitude=151.913773, Latitude=-23.44306 ), proj4string=CRS("+proj=longlat + datum=wgs84"))
Colony_lhi<-SpatialPoints(data.frame(Longitude=159.05991, Latitude=-31.52459 ), proj4string=CRS("+proj=longlat + datum=wgs84"))

DgProj <- CRS("+proj=laea +lon_0=156 +lat_0=-17")

heronProj <- spTransform(Colony_heron, CRS=DgProj)
lhiProj <- spTransform(Colony_lhi, CRS=DgProj)

heronBuffProj <- gBuffer(heronProj, width=80000, quadsegs=50)
lhiBuffProj <- gBuffer(lhiProj, width=80000, quadsegs=50)
#TBuffProj@polygons[[1]]@ID <- as.character(i)
  
heronBuffWgs <- spTransform(heronBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))
lhiBuffWgs <- spTransform(lhiBuffProj, CRS=CRS( "+proj=longlat +ellps=WGS84"))

