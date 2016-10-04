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
daterz=c("mbchla8day", "mhchla8day", "erdVH2chla8day",
                      "agssta8day", "erdMBsstd8day", "jplMURSST41SST",
                       "erdQAstress1daymodStress", "jplMURSST41anom1day")
                      
#"tasshd1day" 2012 max date??

d1<-daterz[1:3]
d2<-daterz[4:6]
d3<-daterz[7:9]

#d1<-dat[dat$Colony=="Heron",]
out<-NULL
for(i in daterz)
  
  tempy<-NULL
  for(j in 1:nrow(dat))
  {
  ext1<-xtracto(xpos=dat$Longitude, ypos=dat$Latitude, tpos=dat$DateAEST,
              dtype=i, xlen=0, ylen=0, verbose=T)
    
  if(class(ext1)=="try-error"){ext1<-xtracto(xpos=dat[j,]$Longitude, ypos=dat[j,]$Latitude, tpos=dat[j,]$DateAEST,
                dtype=i, xlen=0, ylen=0, verbose=T)}
  
  if(class(ext1)=="try-error"){ext1<-errorz}
  
  tempy<-rbind(tempy, ext1)
  }
  
out<-rbind(out, ext1)
print(i)
}

#its well slow.. hmm maybe best using xtracto 3d for monthly mday sst anomoly

xpos <- c(140, 180)
ypos <- c(-5, -43)
tpos <- c('2014-01-01', '2014-05-01') 

test<-xtracto_3D(xpos, ypos, tpos, 'agssta8day', verbose=T)

# trial with rerdap
library(rerddap)
library(raster)
library(ncdf4)
# we're gonna save and export as a .nc file as waaaaaay smaller than csv file size (50 vs 600 Mb)

ed_search(query = 'ssh', which = "grid")$info

out<-info("nrlHycomGLBu008e911S")

(res <- griddap("nrlHycomGLBu008e910S",
                time = c('2014-02-01', '2014-04-08'),
                latitude = c(-10, -42),
                longitude = c(140, 170))) 

(res2 <- griddap("nrlHycomGLBu008e911S",
                time = c('2015-02-01', '2015-04-30'),
                latitude = c(-10, -42),
                longitude = c(140, 170)))

(res3 <- griddap("nrlHycomGLBu008e911S",
                time = c('2016-02-01', '2016-04-18'),
                latitude = c(-10, -42),
                longitude = c(140, 170)))

rezz<-c(res, res2, res3)

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
             vals=d1$surf_el, crs=CRS("+proj=longlat +ellps=WGS84"))
  r1<-flip(r1, direction="y")
  
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
outfile <- paste("ssh_141516.nc")
setwd("~/grive/phd/sourced_data/env_data/erdap_hires")
writeRaster(megastack, outfile, overwrite=TRUE, format="CDF", varname="ssh", varunit="m", 
            longname="ssh -- HYCOM SSH", xname="lon", yname="lat",
            zname="Date", zunit="numeric")

data.nc<- nc_open("ssh_141516.nc")
Zdim = ncvar_get(data.nc,varid="Date")

#print(Zdim)

## Used res for data

#ed_search(query = 'ssh', which = "grid")$info
#info("nrlHycomGLBu008e911S")
#"nrlHycomGLBu008e910S", time = c('2014-02-01', '2014-04-08')
#"nrlHycomGLBu008e911S", time = c('2015-02-01', '2015-04-30')
#"nrlHycomGLBu008e911S",time = c('2016-02-01', '2016-04-18')
             
