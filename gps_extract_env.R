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
                      "tasshd1day", "erdQAstress1daymodStress", "jplMURSST41anom1day",
                      "ETOPO360")

#d1<-dat[dat$Colony=="Heron",]
out<-NULL
for(i in daterz)
  
ext1<-xtracto(xpos=dat$Longitude, ypos=dat$Latitude, tpos=dat$DateAEST,
              dtype=i, xlen=0, ylen=0, verbose=T)
out<-rbind(out, ext1)
print(i)
}

#its well slow.. hmm maybe best using xtracto 3d for monthly mday sst anomoly

xpos <- c(140, 180)
ypos <- c(-5, -43)
tpos <- c('2014-01-01', '2014-05-01') 

test<-xtracto_3D(xpos, ypos, tpos, 'jplMURSST41anommday', verbose=T)



