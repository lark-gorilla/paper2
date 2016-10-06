# 30/09/16 Lancaster, UK
# Extract env data using GPS foraging locations against transiting locations

rm(list=ls())

library(ggplot2)
#devtools::install_github("rmendels/xtractomatic")
library(xtractomatic)

setwd("~/grive/phd/analyses/paper2")

dat<-read.csv("~/grive/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm_attribs.csv", h=T)

#prepdata

dat<-dat[dat$trip_type=="L",]
qplot(Longitude, Latitude, data=dat, colour=factor(Year))

# good no 2014 heron in there.. dont think there are any LT anyway

## Use product code from here:http://coastwatch.pfel.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=griddata
## followed by time unit in dtype parameter, by nchar == 11|12 it allows products not included in original xtracto code to be downloaded

searchData(searchList = list(list("varname", "chl"),
                             list("datasetname", "8day")))

# chl use mbchla8day or mhchla8day or erdVH2chla8day
# sst agssta8day or erdMBsstd8day or jplMURSST41SST
# ssh tasshd1day
# wind erdQAstress1daymodStress
# sst anomaly jplMURSST41anom1day

# should do a extract3d on the blended sst anomaly
"nrlHycomGLBu008e911S", "erdVH3chla8day", "erdMH1chla8day",
  "erdMBsstd8day", "ncdcOisst2Agg",
  "erdQAstress1daymodStress", "jplMURSST41anom1day", "erdAGtanm8day")
                      


# trial with rerdap
library(rerddap)
library(raster)
library(ncdf4)
# we're gonna save and export as a .nc file as waaaaaay smaller than csv file size (50 vs 600 Mb)

ed_search(query = 'erdVH3chla8day', which = "grid")$info

info("ncdcOisst2Agg")

(res <- griddap("ncdcOisst2Agg",
                time = c('2014-02-01', '2014-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170))) 

(res2 <- griddap("ncdcOisst2Agg",
                time = c('2015-02-01', '2015-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170)))

(res3 <- griddap("ncdcOisst2Agg",
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
             vals=d1[,4], crs=CRS("+proj=longlat +ellps=WGS84"))
  
  r1<-flip(r1, direction="y") # turn on for ssh
  
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
outfile <- paste("sstA_ncdcOisst2Agg.nc")
setwd("~/grive/phd/sourced_data/env_data/erdap_hires")
writeRaster(megastack, outfile, overwrite=TRUE, format="CDF", varname="sstA", varunit="deg C", 
            longname="SST anomoly -- ncdcOisst2Agg", xname="lon", yname="lat",
            zname="Date", zunit="numeric")

data.nc<- nc_open("chla_erdVH3chla8day.nc")
Zdim = ncvar_get(data.nc,varid="Date")

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
