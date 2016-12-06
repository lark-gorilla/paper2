  # 13/09/16 Lancaster, UK
# Clean and compile all GPS tracking data for use in multiple PhD chapters
# Some things to remember:
# Date, Time and DateTime are the unresampled Date and Time values (AEST),
# Tracktime, DateAEST and TimeAESTare resamapled, These should 
# be used for analyses, the earlier ones are for reference.

rm(list=ls())

library(maps)
library(sp)
library(geosphere)

setwd("~/grive/phd/analyses/tracking_data_pot")

# read in heron 2014 resamp later

heron_2015<-read.csv("~/grive/phd/sourced_data/Heron/Heron_2015_GPS/Heron_2015_GPS_points.csv", h=T)

lhi_2014<-read.csv("~/grive/phd/fieldwork/LHI_Mar_2014/data/tracking_results/compiled_tracking_data.csv", h=T)

# LHI2014 suffers from UTC time issue (all points in UTC), need to
# correct before proceeding
lhi_2014$DateTime<-paste(lhi_2014$Date, lhi_2014$Time, sep="")
lhi_2014$DateTime <- as.POSIXct(strptime(lhi_2014$DateTime, "%Y/%m/%d %H:%M:%S"), "GMT")
lhi_2014$TrackTime <- as.double(lhi_2014$DateTime)
# bit below advances time to UTC +11
lhi_2014$TrackTime<-lhi_2014$TrackTime+(11*60)*60
lhi_2014$Date <- as.Date(as.POSIXlt(lhi_2014$TrackTime, origin="1970-01-01", "GMT"))
lhi_2014$Time <- format((as.POSIXlt(lhi_2014$TrackTime, origin="1970-01-01", "GMT")), "%H:%M:%S")

lhi_2015<-read.csv("~/grive/phd/fieldwork/LHI_Feb_2015/data/tracking_results/compiled_tracking_data.csv", h=T)
lhi_2016<-read.csv("~/grive/phd/fieldwork/LHI_Feb_2016/data/tracking_data/compiled_tracking_data.csv", h=T)
                            
#makes sure we're using the AEST datetime NOT the UTC one!
heron_2015$Date<-paste(substr(heron_2015$DateTime_AEST, 1, 4),"/",
                       substr(heron_2015$DateTime_AEST, 6, 7),"/",
                       substr(heron_2015$DateTime_AEST, 9, 10), sep="")
heron_2015$Time<-substr(heron_2015$DateTime_AEST, 12, 19)


dat<-rbind(data.frame( TrackID=heron_2015$TrackID, Date=heron_2015$Date, Time=heron_2015$Time, 
                      Latitude=heron_2015$Latitude, Longitude=heron_2015$Longitude, Colony="Heron"), 
           data.frame(TrackID=lhi_2014$TrackID, Date=lhi_2014$Date, Time=lhi_2014$Time,
                      Latitude=lhi_2014$Latitude, Longitude=lhi_2014$Longitude,  Colony="LHI"),
           data.frame(TrackID=lhi_2015$TrackID, Date=lhi_2015$Date, Time=lhi_2015$Time,
                      Latitude=lhi_2015$Latitude, Longitude=lhi_2015$Longitude, Colony="LHI"),
           data.frame(TrackID=lhi_2016$TrackID, Date=lhi_2016$Date, Time=lhi_2016$Time,
                      Latitude=lhi_2016$Latitude, Longitude=lhi_2016$Longitude, Colony="LHI"))

dat$DateTime<-paste(dat$Date, dat$Time, sep="")
dat$DateTime <- as.POSIXct(strptime(dat$DateTime, "%Y/%m/%d %H:%M:%S"), "GMT")

dat$TrackTime <- as.double(dat$DateTime)
head(dat)


source("~/research/seabird_analyses/Heron_Island/heron_analyses_r_GPS/Velocity_Functions.r")

tripz<-unique(dat$TrackID)

###@@@ SPEED FILTER @@@###

datout <- NULL
for(i in tripz)
{
  Temp <-  dat[dat$TrackID == i,]
  plot(Latitude~Longitude, data=Temp, asp=1, cex=0.5, col=4, main=Temp$Bird_ID[1])
  map("world", add=T)
  
  Temp$Vel <- 0
  for(j in 1:nrow(Temp))
  {
    Temp[j,]$Vel <- backforVel(point=j, trip=Temp, n=4, alt=FALSE, filter=FALSE)
    plot(Latitude~Longitude, data=Temp, asp=1, cex=0.5, col=4)
  }
  points(Latitude~Longitude, data=Temp[Temp$Vel > 75,], asp=1, cex=0.5, col=2)
  datout <- rbind(datout, Temp)
  #readline("OK")
}

###@@@ removing points going with speeds over 75 km/h @@@###
datout<-datout[datout$Vel<75,]
datout<-na.omit(datout)

####@@@@ RESAMPLE @@@@####

source("~/research/seabird_analyses/Heron_Island/heron_analyses_r_GPS/Resample.r")

## set up to not interpolate between points > 1 hr apart

results<-NULL
for(i in tripz)
{
  Track<-datout[datout$TrackID == i,]
  
  #Track<-Track[-which(duplicated(Track$TrackTime)==TRUE),]
  #print(paste(length(which(duplicated(Track$TrackTime)==TRUE)), "duplicates removed", sep=" "))
  
  resample_output<-resample(Track, timeStep=0.1666667)  ## timeStep set for 10 min intervals
  
  #readline("ok")
  results<-rbind(results,resample_output)
  print(i)
}

plot(Latitude~Longitude, results, pch=16, cex=0.4, col=Bird_ID)
map("worldHires", add=T, col=3)

# IMPORTANT! we call these AEST as they are, but in the code we specify GMT to keep R happy
results$DateAEST <- as.Date(as.POSIXlt(results$TrackTime, origin="1970-01-01", "GMT"))
results$TimeAEST <- format((as.POSIXlt(results$TrackTime, origin="1970-01-01", "GMT")), "%H:%M:%S")

results$DateTime2 <- paste(results$DateAEST, results$TimeAEST, sep= " ")
results$DateTime2 <- as.POSIXct(strptime(results$DateTime2, "%Y-%m-%d %H:%M:%S"), "GMT")
results$TrackTime2 <- as.double(results$DateTime2)

write.csv(results, "GPS_141516_clean_resamp.csv", quote=F, row.names=F)

results<-read.csv("GPS_141516_clean_resamp.csv", h=T)
results<-results[results$Longitude>147,] # kill erroneous points
results$DateTime2<-NULL
results$TrackTime2<-NULL

heron_2014<-read.csv("~/grive/phd/sourced_data/Heron/GPS_ST_2014_clean_resamp.csv", h=T)
#this data is already speed filtered and interpolated to 10 min interval 

heron_2014<-data.frame(Date=heron_2014$Date, Time=heron_2014$time,
                       Latitude=heron_2014$Latitude, Longitude=heron_2014$Longitude, 
                       TrackID=paste(heron_2014$Nest_ID, heron_2014$Bird_ID, sep="_"), Colony="Heron",
                       DateTime=heron_2014$DateTime, TrackTime=heron_2014$TrackTime, 
                       Vel=heron_2014$Vel, DateAEST=heron_2014$DateGMT, TimeAEST=heron_2014$TimeGMT)

results<-rbind(results, heron_2014)
write.csv(results, "GPS_141516_clean_resamp.csv", quote=F, row.names=F)


### Colony Lookup for tripsplit ###

Colony_heron<-data.frame(Latitude=-23.44306, Longitude=151.913773)
Colony_lhi<-data.frame(Latitude=-31.52459, Longitude=159.05991)


#### Trip split ####
source("~/grive/phd/scripts/MIBA_scripts_revised/TripSplit_revised_MMLINUX.r")
library(maptools)

results$ID<-results$TrackID

for(i in unique(results[results$Colony=="Heron",]$ID))
{
  Temp <- results[results$TrackID==i,]
  
  Trip <- tripSplit(Track=Temp, Colony=Colony_heron, InnerBuff=1, ReturnBuff=20, Duration=10, plotit=T)
  
  if(which(unique(results[results$Colony=="Heron",]$ID)==i) == 1) {Trips <- Trip} else
    Trips <- spRbind(Trips,Trip)
  
  #readline("")
}


plot(Trips)
plot(Trips[Trips$trip_id==-1,], col=2, add=T)
Trips<-Trips[Trips$trip_id!=-1,]

Trips@data$PointId<-0
for(i in unique(Trips$trip_id))
{
  Trips@data[Trips$trip_id==i,]$PointId<-seq(1:nrow(Trips@data[Trips$trip_id==i,]))
}

source("~/research/phd/scripts/MIBA_scripts_revised/TripSplit_revised_MMLINUX.r")

dat_tripsplit<-Trips@data # store heron tripsplit

for(i in unique(results[results$Colony=="LHI",]$ID))
{
  Temp <- results[results$TrackID==i,]
  
  Trip <- tripSplit(Track=Temp, Colony=Colony_lhi, InnerBuff=1, ReturnBuff=20, Duration=10, plotit=T)
  
  if(which(unique(results[results$Colony=="LHI",]$ID)==i) == 1) {Trips <- Trip} else
    Trips <- spRbind(Trips,Trip)
  
  #readline("")
}


plot(Trips)
plot(Trips[Trips$trip_id==-1,], col=2, add=T)
## removes -1s
Trips<-Trips[Trips$trip_id!=-1,]

Trips@data$PointId<-0
for(i in unique(Trips$trip_id))
{
  Trips@data[Trips$trip_id==i,]$PointId<-seq(1:nrow(Trips@data[Trips$trip_id==i,]))
}

dat_tripsplit<-rbind(dat_tripsplit, Trips@data)

# write out tripsplit data

write.csv(dat_tripsplit, "GPS_141516_clean_resamp_tripsplit.csv", quote=F, row.names=F)

## HMM trial using moveHMM package
library(moveHMM)

dat<-read.csv("GPS_2014_15_clean_resamp_tripsplit.csv", h=T)

tempdata<-prepData(data.frame(x=dat$Longitude, y=dat$Latitude,
                              ID=dat$TrackID), type="LL")
# could use Colony and ColDist as covariates
plot(tempdata) # this gives idea of starting parameter values
plot(tempdata, compact=T) #better plot

mu0 <- c(0.5,2,5)
sigma0 <- c(0.5,2,5)
zeromass0 <- c(0.5,0.1,0.1)
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0,pi,0)
kappa0 <- c(2,1,2)
anglePar0 <- c(angleMean0,kappa0)

hmm_1<- fitHMM(data=tempdata,nbStates=3,stepPar0=stepPar0,
                anglePar0=anglePar0)

states<-viterbi(hmm_1)
dat$HMMstates<-states
write.csv(dat, "GPS_141516_clean_resamp_tripsplit_hmm.csv", quote=F, row.names=F)

# fix to correct UTC error times in lhi 2014 data
dat<-read.csv("GPS_141516_clean_resamp_tripsplit_hmm.csv", h=T)
dat$Year<-substr(dat$Date, 1,4)

d1<-dat[which(dat$Colony=="LHI" & dat$Year==2014),]
d2<-dat[-which(dat$Colony=="LHI" & dat$Year==2014),]

d1$TrackTime<-d1$TrackTime+(11*60)*60

d1$DateAEST<-as.character(d1$DateAEST)

d1$DateAEST<-as.Date(as.POSIXlt(d1$TrackTime,origin="1970-01-01", "GMT"))
d1$TimeAEST <-format((as.POSIXlt(d1$TrackTime,origin="1970-01-01", "GMT")), "%H:%M:%S")

dat<-rbind(d1, d2)

write.csv(dat, "GPS_141516_clean_resamp_tripsplit_hmm.csv", quote=F, row.names=F)

### Now we make the trip summary statistics

### SUMMARISE MAX DIST FROM COLONY AND TRIP TRAVELLING TIME FOR EACH TRIP

trip_distances<-data.frame(trip=unique(dat$trip_id), max_dist=0, time=0, total_dist=0, fin_dist=0,
                           Returns="na",trip_type="na", Year="na", Colony="na",
                           perc_forage=0, perc_commute=0, perc_sit=0)
                           #,perc_forage_day=0, perc_sit_day=0)  	### create data frame for each trip

trip_distances$Returns<-as.character(trip_distances$Returns)
trip_distances$trip_type<-as.character(trip_distances$trip_type)
trip_distances$Year<-as.character(trip_distances$Year)
trip_distances$Colony<-as.character(trip_distances$Colony)

library(geosphere)
for (i in trip_distances$trip)
  {
  trip_distances[trip_distances$trip==i,2]<-max(dat[dat$trip_id==i,]$ColDist)/1000
  trip_distances[trip_distances$trip==i,3]<-(max(dat[dat$trip_id==i,]$TrackTime)-min(dat[dat$trip_id==i,]$TrackTime))/3600
  trip_distances[trip_distances$trip==i,5]<-max(dat[dat$trip_id==i,]$ColDist)/1000
  ## Calculate distances from one point to the next and total trip distance
  x=dat[dat$trip_id==i,]
  x$Dist=0
  x$Dist[1]<-x$ColDist[1]/1000				### distance to first point is assumed a straight line from the nest/colony
  for (p in 2:dim(x)[1])
    {
    p1<-c(x$Longitude[p-1],x$Latitude[p-1])
    p2<-c(x$Longitude[p],x$Latitude[p])
    #x$Dist[p]<-pointDistance(p1,p2, lonlat=T, allpairs=FALSE)/1000			### no longer works in geosphere
    x$Dist[p]<-distMeeus(p1,p2)/1000						### great circle distance according to Meeus, converted to km
    }
  trip_distances[trip_distances$trip==i,4]<-sum(x$Dist)+(x$ColDist[p]/1000)	## total trip distance is the sum of all steps plus the dist from the nest of the last location - for non return trips this will be an underestimate
  trip_distances[trip_distances$trip==i,6]<-as.character(unique(dat[dat$trip_id==i,]$Returns))
  trip_distances[trip_distances$trip==i,7]<-"L"
  trip_distances[trip_distances$trip==i,8]<-as.character(substr(dat[dat$trip_id==i,]$DateAEST[1],1,4))
  trip_distances[trip_distances$trip==i,9]<-as.character(unique(dat[dat$trip_id==i,]$Colony))
  trip_distances[trip_distances$trip==i,10]<-round((nrow(dat[dat$trip_id==i &dat$HMMstates==2,])/
                                                      nrow(dat[dat$trip_id==i, ]))*100)
  trip_distances[trip_distances$trip==i,11]<-round((nrow(dat[dat$trip_id==i &dat$HMMstates==3,])/
                                                      nrow(dat[dat$trip_id==i, ]))*100)
  trip_distances[trip_distances$trip==i,12]<-round((nrow(dat[dat$trip_id==i &dat$HMMstates==1,])/
                                                      nrow(dat[dat$trip_id==i, ]))*100)

  #trip_distances[trip_distances$trip==i,13]<-round((nrow(dat[dat$trip_id==i &dat$HMMstates==2 & dat$DayNight=="Day",])/
                                                      #nrow(dat[dat$trip_id==i &dat$HMMstates==2,]))*100)
  #trip_distances[trip_distances$trip==i,14]<-round((nrow(dat[dat$trip_id==i &dat$HMMstates==1 & dat$DayNight=="Day",])/
                                                      #nrow(dat[dat$trip_id==i &dat$HMMstates==1,]))*100)
  }

trip_distances$days<-ceiling(trip_distances$time/24)

#correctly specify short and long trips
trip_distances[trip_distances$time<24*4,]$trip_type<-"S"
trip_distances[trip_distances$trip_type=="S" &
                 trip_distances$Returns=="N" &
                 trip_distances$fin_dist>202,]$trip_type<-"L"
# this is a kind of best guess but need to validate these 'long' trips

#write out trip distances
trip_distances[trip_distances$Returns=="",]$Returns<-"Y"
write.csv(trip_distances, "GPS_141516_trip_summary.csv", quote=F, row.names=F)


### Add some additional trip descriptor columns and write out
dat$days<-trip_distances[match(dat$trip_id,trip_distances$trip),]$days
dat$trip_type<-trip_distances[match(dat$trip_id,trip_distances$trip),]$trip_type
# fix for some cheeky little LTs that are actually STs
dat[dat$trip_id %in% c("25_02_15_06LW3", "28_03_15_04LW4",
                       "16_02_16_41LW3", "18_02_16_22LW3"),]$trip_type<-"S"


dat$Month<-substr(dat$DateAEST, 6,7)
dat$Year<-substr(dat$DateAEST, 1,4)
dat$Hr<-substr(dat$TimeAEST, 1,2)

#the below code accurately assigns day and night to each datapoint 
library(RAtmosphere)
# this code is ripped from https://smathermather.wordpress.com/2014/08/05/using-spatial-data-in-r-to-estimate-home-ranges-guest-blog-post/
  #if time is greater than the sunrise time or less than the sunset time, apply DAY
  suntime <- suncalc(as.numeric(as.Date(dat$DateAEST, format="%Y-%m-%d") - as.Date(paste(dat$Year,"-01-01", sep=""))), Lat=dat$Latitude, Long=dat$Longitude, UTC=FALSE)
  coytime <- sapply(strsplit(as.character(dat$TimeAEST),":"),
                    function(x) {
                      x <- as.numeric(x)
                      x[1]+x[2]/60+x[3]/3600
                    })
  
  dat$DayNight<-ifelse(coytime > suntime$sunrise & coytime < suntime$sunset, "DAY", "NIGHT")

write.csv(dat, "GPS_141516_clean_resamp_tripsplit_hmm_attribs.csv", quote=F, row.names=F)

# add in extra trip descriptor fields
dat$DateTimeAEST<-paste(dat$DateAEST, dat$TimeAEST, sep="")
dat$DateTimeAEST <- as.POSIXct(strptime(dat$DateTimeAEST, "%Y-%m-%d %H:%M:%S"), "GMT")

trip_distances$stDateTime<-min(dat$DateTimeAEST)
trip_distances$edDateTime<-min(dat$DateTimeAEST)
for (i in trip_distances$trip)
{
  trip_distances[trip_distances$trip==i,]$stDateTime<-
    min(dat[dat$trip_id==i,]$DateTimeAEST)
  trip_distances[trip_distances$trip==i,]$edDateTime<-
    max(dat[dat$trip_id==i,]$DateTimeAEST)
print(i)
}

write.csv(trip_distances, "GPS_141516_trip_summary.csv", quote=F, row.names=F)
# read back in now with return dates added
trip_distances<-read.csv("GPS_141516_trip_summary.csv", h=T) 

stDateTime<-as.POSIXlt(trip_distances$stDateTime, format = "%Y-%m-%d %H:%M:%S")

DateTime_AR<-as.POSIXlt(trip_distances$Actual_return, format = "%d/%m/%Y %H:%M")

trip_distances$Actual_time<-difftime(DateTime_AR,stDateTime, units="hours")

write.csv(trip_distances, "GPS_141516_trip_summary.csv", quote=F, row.names=F)
# re



