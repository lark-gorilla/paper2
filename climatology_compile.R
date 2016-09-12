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