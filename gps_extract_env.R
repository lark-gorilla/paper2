# 30/09/16 Lancaster, UK
# Extract env data using GPS foraging locations against transiting locations

rm(list=ls())
# rerdap
library(rerddap)
library(raster)
library(ncdf4)

# We use rerddap package to make direct calls to the erddap server and get 
# gridded data via the OPeNDAP hyperslab protocol??

#product code from here:
#https://coastwatch.pfeg.noaa.gov/erddap/griddap/index.html?page=1&itemsPerPage=1000

# chl use erdVH3chla8day or erdMH1chla8day 
# sst erdAGssta8day or ncdcOisst2Agg
# ssh nrlHycomGLBu008e911S
# wind stress and upwelling erdQMstress3day
# sst anomaly ncdcOisst2Agg or erdAGtanm8day

# should do a extract3d on the blended sst anomaly
"nrlHycomGLBu008e911S", "erdVH3chla8day", "erdMH1chla8day",
  "erdAGssta8day", "ncdcOisst2Agg", "erdAGtanm8day"
  "erdQMstress8day" 
                      
# we're gonna save and export as a .nc file as waaaaaay smaller than csv file size (50 vs 600 Mb)

ed_search(query = 'erdQMstress1day', which = "grid")$info

info("erdQMstress1day")

(res <- griddap("erdQMstress8day",
                time = c('2014-02-01', '2014-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170))) 

(res2 <- griddap("erdQMstress8day",
                time = c('2015-02-01', '2015-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170)))

(res3 <- griddap("erdQMstress8day",
                time = c('2016-02-01', '2016-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170)))

for( j in 1:3)
{
k<-NULL  
if(j==1){k<-res}
  if(j==2){k<-res2}
  if(j==3){k<-res3}
  
for( i in unique(k$data$time))
{
  d1<-k$data[k$data$time==i,]
  r1<-raster(xmn=min(d1$lon), xmx=max(d1$lon), ymn=min(d1$lat), ymx=max(d1$lat),
             nrows=length(unique(d1$lat)), ncols=length(unique(d1$lon)), 
             vals=d1[,8], crs=CRS("+proj=longlat +ellps=WGS84"))
            # the correct column needs to be set for vals, normally 4
  r1<-flip(r1, direction="y") # turn off for chla
  
  if(which(i==unique(k$data$time))==1){st1<-r1}else{
    st1<-stack(st1, r1)}
 print(i)
}

if(j==1){megastack<-st1}else{
  megastack<-stack(megastack, st1)}
}

megastack<-setZ(megastack, 
                c(paste(substr(unique(res$data$time), 1,4),
                     substr(unique(res$data$time), 6,7),
                     substr(unique(res$data$time), 9,10),sep=""),
                  paste(substr(unique(res2$data$time), 1,4),
                        substr(unique(res2$data$time), 6,7),
                        substr(unique(res2$data$time), 9,10),sep=""),
                  paste(substr(unique(res3$data$time), 1,4),
                        substr(unique(res3$data$time), 6,7),
                        substr(unique(res3$data$time), 9,10),sep="")))


# Save the raster file as a netCDF
outfile <- paste("ekmU_erdQMstress8day.nc")
setwd("~/grive/phd/sourced_data/env_data/erdap_hires")
writeRaster(megastack, outfile, overwrite=TRUE, format="CDF", varname="ekeman upwelling", varunit="m s-1", 
            longname="Wind -- erdQMstress8day", xname="lon", yname="lat",
            zname="Date", zunit="numeric")

data.nc<- nc_open("sst_ncdcOisst2Agg.nc")
Zdim = ncvar_get(data.nc,varid="Date")

#make visualisation

library(rasterVis)
library(animation)
library(gridExtra)
setwd("~/grive/phd/sourced_data/env_data/erdap_hires")

data.nc<- nc_open("sst_ncdcOisst2Agg.nc")
Zdim = ncvar_get(data.nc,varid="Date")
r2014<-grep("2014", Zdim)
r2015<-grep("2015", Zdim)
r2016<-grep("2016", Zdim)

lhi16<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2016.shp", layer="LTLHIGPS2016")
lhi15<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2015.shp", layer="LTLHIGPS2015")
lhi14<-readOGR(dsn="/home/mark/grive/phd/analyses/paper2/spatial/LTLHIGPS2014.shp", layer="LTLHIGPS2014")

saveGIF({
for(i in 1:89){
r1<-raster("sst_ncdcOisst2Agg.nc", band=r2014[i])  
r2<-raster("sst_ncdcOisst2Agg.nc", band=r2015[i])
r3<-raster("sst_ncdcOisst2Agg.nc", band=r2016[i])
  
col.l <- colorRampPalette(c('black','blue','cyan','yellow', 'orange', 'red'))(20)
p1<-levelplot(r1, contour=T, margin=F,scales=list(draw=FALSE) , at=13:33, col.regions=col.l, ylab= NULL, xlab= NULL, colorkey=T,main=paste(r1@z[[1]]))+layer(sp.polygons(lhi14[c(1,3),], col=c("dark gray","red")))
p2<-levelplot(r2, contour=T, margin=F,scales=list(draw=FALSE) ,ylab= NULL,at=13:33,col.regions=col.l, xlab= NULL,colorkey=T, main=paste(r2@z[[1]]))+layer(sp.polygons(lhi15[c(1,3),], col=c("dark gray","red")))
p3<-levelplot(r3, contour=T, margin=F,scales=list(draw=FALSE) ,at=13:33,col.regions=col.l,ylab= NULL, xlab= NULL,colorkey=T, main=paste(r3@z[[1]]))+layer(sp.polygons(lhi16[c(1,3),], col=c("dark gray","red")))
grid.arrange(p1,p2,p3, ncol=3)}

}, interval=1.5, "sst_141516V4.gif", ani.width=1500, ani.height=500)

#print(Zdim)

## Used res for data

#ed_search(query = 'ssh', which = "grid")$info
#info("nrlHycomGLBu008e911S")
#"nrlHycomGLBu008e910S", time = c('2014-02-01', '2014-04-08')
#"nrlHycomGLBu008e911S", time = c('2015-02-01', '2015-04-30')
#"nrlHycomGLBu008e911S",time = c('2016-02-01', '2016-04-18')
 
#"erdMH1chla8day", time = c('2014-02-01', '2014-04-30'),
#"erdMH1chla8day",time = c('2015-02-01', '2015-04-30'),
#"erdMH1chla8day",time = c('2016-02-01', '2016-04-30'),
                
#erdVH3chla8day
#"erdVH3chla8day", time = c('2014-02-01', '2014-04-30'),
#"erdVH3chla8day",time = c('2015-02-01', '2015-04-30'),
#"erdVH3chla8day",time = c('2016-02-01', '2016-04-30'), # end date early

#ncdcOisst2Agg # we used sst and anom products from this variable
#"ncdcOisst2Agg", time = c('2014-02-01', '2014-04-30'),
#"ncdcOisst2Agg",time = c('2015-02-01', '2015-04-30'),
#"ncdcOisst2Agg",time = c('2016-02-01', '2016-04-30'),

#erdAGssta8day
#"erdAGssta8day", time = c('2014-02-01', '2014-04-30'),
#"erdAGssta8day",time = c('2015-02-01', '2015-04-30'),
#"erdAGssta8day",time = c('2016-02-01', '2016-04-17'), end of timelimited

#erdQMstress1day # we used wind stress modulus and ekman upwelling products from this variable
#"erdQMstress1day", time = c('2014-02-01', '2014-04-30'),
#"erdQMstress1day",time = c('2015-02-01', '2015-04-30'),
#"erdQMstress1day",time = c('2016-02-01', '2016-04-30'),


# now to extract the data!
rm(list=ls())
library(raster)
library(ncdf4)
library(ggplot2)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("~/grive/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm_attribs.csv", h=T)

#prepdata

dat<-dat[dat$trip_type=="L",]
qplot(Longitude, Latitude, data=dat, colour=factor(Year))
# good no 2014 heron in there.. dont think there are any LT anyway

oceo<-list.files("~/grive/phd/sourced_data/env_data/erdap_hires")

d_chlMH1<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/chla_erdMH1chla8day.nc")
d_chlVH3<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/chla_erdVH3chla8day.nc")
d_ekmU<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/ekmU_erdQMstress8day.nc")
d_modW<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/modW_erdQMstress8day.nc")
d_sstAG<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/sst_erdAGssta8day.nc")
d_sstOi<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/sst_ncdcOisst2Agg.nc")
d_AsstAG<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/sstA_erdAGtanm8day.nc")
d_AsstOi<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/sstA_ncdcOisst2Agg.nc")
d_sshHy<- nc_open("~/grive/phd/sourced_data/env_data/erdap_hires/ssh_nrlHycomGLBu008e911S.nc")

d_chlMH1_D<-ncvar_get(d_chlMH1,varid="Date") 
d_chlVH3_D<- ncvar_get(d_chlVH3,varid="Date")
d_ekmU_D<-ncvar_get(d_ekmU,varid="Date")
d_modW_D<- ncvar_get(d_modW,varid="Date")
d_sstAG_D<- ncvar_get(d_sstAG,varid="Date")
d_sstOi_D<- ncvar_get(d_sstOi,varid="Date")
d_AsstAG_D<- ncvar_get(d_AsstAG,varid="Date")
d_AsstOi_D<- ncvar_get(d_AsstOi,varid="Date")
d_sshHy_D<- ncvar_get(d_sshHy,varid="Date")

dat$chlMH1<-0
dat$chlVH3<-0
dat$ekmU<-0
dat$modW<-0
dat$sstAG<-0
dat$sstOi<-0
dat$AsstAG<-0
dat$AsstOi<-0
dat$sshHy<-0

setwd("~/grive/phd/sourced_data/env_data/erdap_hires")

for( i in unique(dat$DateAEST))
 {
 iTT<-as.double(as.POSIXct(strptime(i, "%Y-%m-%d"), "GMT"))
 
 #sshHy
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_sshHy_D, "%Y%m%d"), "GMT")))^2))==
         sqrt((iTT-as.double(as.POSIXct(strptime(d_sshHy_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$sshHy<-extract(raster("ssh_nrlHycomGLBu008e911S.nc", band=k1[1]),
                                      data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                 dat[dat$DateAEST==i,]$Latitude))
 #chlMH1
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_chlMH1_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_chlMH1_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$chlMH1<-extract(raster("chla_erdMH1chla8day.nc", band=k1[1]),
                                      data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                 dat[dat$DateAEST==i,]$Latitude))
 #chlVH3
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_chlVH3_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_chlVH3_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$chlVH3<-extract(raster("chla_erdVH3chla8day.nc", band=k1[1]),
                                       data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                  dat[dat$DateAEST==i,]$Latitude))
 #ekmU
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_ekmU_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_ekmU_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$ekmU<-extract(raster("ekmU_erdQMstress8day.nc", band=k1[1]),
                                       data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                  dat[dat$DateAEST==i,]$Latitude))
 #modW
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_modW_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_modW_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$modW<-extract(raster("modW_erdQMstress8day.nc", band=k1[1]),
                                     data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                dat[dat$DateAEST==i,]$Latitude))
 #sstAG
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_sstAG_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_sstAG_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$sstAG<-extract(raster("sst_erdAGssta8day.nc", band=k1[1]),
                                     data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                dat[dat$DateAEST==i,]$Latitude))
 #sstOi
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_sstOi_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_sstOi_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$sstOi<-extract(raster("sst_ncdcOisst2Agg.nc", band=k1[1]),
                                      data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                 dat[dat$DateAEST==i,]$Latitude))
 #AsstAG
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_AsstAG_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_AsstAG_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$AsstAG<-extract(raster("sstA_erdAGtanm8day.nc", band=k1[1]),
                                      data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                 dat[dat$DateAEST==i,]$Latitude))
 #AsstOi
 k1<-which(min(sqrt((iTT-as.double(as.POSIXct(strptime(d_AsstOi_D, "%Y%m%d"), "GMT")))^2))==
             sqrt((iTT-as.double(as.POSIXct(strptime(d_AsstOi_D, "%Y%m%d"), "GMT")))^2))
 
 dat[dat$DateAEST==i,]$AsstOi<-extract(raster("sstA_ncdcOisst2Agg.nc", band=k1[1]),
                                       data.frame(dat[dat$DateAEST==i,]$Longitude,
                                                  dat[dat$DateAEST==i,]$Latitude))
 print(i)
 }
 
 write.csv(dat, "~/grive/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm_oceano_attribs.csv", quote=F, row.names=F)
 

