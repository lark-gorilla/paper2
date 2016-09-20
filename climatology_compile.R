# 11/09/16 Lancaster, UK
#read in, select and compile hdf data for kernal-based extraction

rm(list=ls())
setwd("~/grive/phd/sourced_data/env_data/climatologies")
library("curl")
#Download data from bloomwatch 360 using get query
# REMEMBER you can grab the get query link by right clicking
# once parameters have been set the .nc download option

#http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=TNCwwww&timePeriod=1month&centeredTime=1995-12-30T00:00:00&minLon=0&maxLon=360&minLat=-90&maxLat=90&fileType=.nc
#http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=CPCsstn&timePeriod=monthly&centeredTime=0001-12-16T12:00:00&minLon=0&maxLon=360&minLat=-90&maxLat=90&fileType=.nc
#http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jThe Changjiang River discharge affects the distribution of foraging seabirdssp?get=gridData&dataSet=CTCsshd&timePeriod=monthly&centeredTime=0001-12-16T12:00:00&minLon=0&maxLon=360&minLat=-90&maxLat=90&fileType=.nc
#http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=CQCwekm&timePeriod=monthly&centeredTime=0001-12-16T12:00:00&minLon=0&maxLon=360&minLat=-90&maxLat=90&fileType=.nc
#http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=CPMbfp1&timePeriod=monthly&centeredTime=0001-12-16T12:00:00&minLon=0&maxLon=360&minLat=-90&maxLat=90&fileType=.nc

dset=c("CPCsstn","CTCsshd","CQCwekm", "CPMbfp1", "TNCwwww")

for(i in dset)
  {
  if(i=="TNCwwww")
    {
    for(j in c("01-29","03-01","03-31", "04-30", "05-30","06-30", "07-30",
               "08-30", "09-30", "10-30", "11-30", "12-30")){
    curl_download(url=paste("http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=TNCwwww&timePeriod=1month&centeredTime=1995-",
                            j, "T00:00:00&minLon=100&maxLon=180&minLat=-45&maxLat=45&fileType=.nc", sep=""), destfile=paste(i, "1995-", j, ".nc", sep=""))
    print(paste(i, j))}
    
    }else{
      for(j in c("01","02","03", "04", "05","06", "07",
                  "08", "09", "10", "11", "12")){
    curl_download(url=paste("http://coastwatch.pfeg.noaa.gov/coastwatch/CWBrowserWW360.jsp?get=gridData&dataSet=",
    i,"&timePeriod=monthly&centeredTime=~0001-", j, "-16T12:00:00&minLon=100&maxLon=180&minLat=-45&maxLat=45&fileType=.nc",sep=""), destfile=paste(i, "2016-", j, "-01", ".nc", sep=""))
        print(paste(i, j))} # the curl above has a ~ before centered time to approx next best time, need to check it got it right
    }
  }  
 
library(ncdf4)
library(raster)
ncz<-list.files()

nc_open(ncz[1])

ras<-raster(ncz[1])

# now get thermocline data
nc1<-nc_open("~/grive/phd/sourced_data/env_data/thermocline/ttd_DTm02_c1m_reg2.0.nc")
print(nc1)

ttd=ncvar_get(nc1, varid="ttd") # extracts top of thermocline depth (estimated by kriging of ttd_smth) data
## from .nc file note gives XYZ matrix of values for long, lats and months [1:180, 1:90, 12]

for(i in 1:12){
r1<-raster(nrows=90, ncols=180, resolution=2,
           vals=as.vector(ttd[1:180, 1:90, i])) # 3 is cos we want March
r1@extent<-extent(0,360,-90,90) # global coords
r1<-flip(r1,2)
r1<-rotate(r1)
#plot(r1)
#plot(log(r1))
r1[r1==1.000000e+09]<-NA
#plot(r1)
writeRaster(r1, paste("~/grive/phd/sourced_data/env_data/climatologies/thermocline_",
              i, ".tif", sep=""))
}

### finally we better make a distance to seamount layer
bathy<-raster("~/grive/phd/sourced_data/env_data/phd_bathy/GRIDONE_2D_100.0_-45.0_180.0_40.0.nc")
smounts<-read.table("~/grive/phd/sourced_data/env_data/seamounts/seamounts_KWSMTSv01.txt",
           skip=16, header=T, sep = "\t", comment.char=">")
plot(Latitude~Longitude, smounts, pch=7, cex=0.5)

sm_ras_wgs<-rasterize(data.frame(smounts$Longitude,
                      smounts$Latitude),bathy, field=1, background=NA) 

sm_ras_merc<-projectRaster(sm_ras_wgs, crs="+proj=merc +lon_0=140 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#project to mercator with centre at 140 

writeRaster(sm_ras_merc, "~/grive/phd/sourced_data/env_data/seamounts/seamounts_ras_merc.tif", overwrite=T) 
## distance (even when projected) in R is pretty slow, if you write to QGIS to use GDAL = v. fast :)

#Now read in QGIS product and reproject
dsmount_merc<-raster("~/grive/phd/sourced_data/env_data/seamounts/d_seamounts_merc.tif")

dsmount_wgs<-projectRaster(dsmount_merc, crs="+proj=longlat +ellps=WGS84")
#convert to km
dsmount_wgs<-dsmount_wgs/1000

writeRaster(dsmount_wgs, "~/grive/phd/sourced_data/env_data/seamounts/d_seamounts_wgs.tif", overwrite=T) 


