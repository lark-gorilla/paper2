rm(list=ls())

library(maps)
library(sp)
library(geosphere)

setwd("~/grive/phd/analyses/tracking_data_pot")

# read in heron 2014 resamp later

heron_2015<-read.csv("~/grive/phd/sourced_data/Heron/Heron_2015_GPS/Heron_2015_GPS_points.csv", h=T)

lhi_2014<-read.csv("~/grive/phd/fieldwork/LHI_Mar_2014/data/tracking_results/compiled_tracking_data.csv", h=T)

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


results$DateGMT <- as.Date(as.POSIXlt(results$TrackTime, origin="1970-01-01", "GMT"))
results$TimeGMT <- format((as.POSIXlt(results$TrackTime, origin="1970-01-01", "GMT")), "%H:%M:%S")

results$DateTime2 <- paste(results$DateGMT, results$TimeGMT, sep= " ")
results$DateTime2 <- as.POSIXct(strptime(results$DateTime2, "%Y-%m-%d %H:%M:%S"), "GMT")
results$TrackTime2 <- as.double(results$DateTime2)

write.csv(results, "GPS_2014_15_clean_resamp.csv", quote=F, row.names=F)

results<-read.csv("GPS_2014_15_clean_resamp.csv", h=T)
results<-results[results$Longitude>147,] # kill erroneous points
results$DateTime2<-NULL
results$TrackTime2<-NULL

heron_2014<-read.csv("~/grive/phd/sourced_data/Heron/GPS_ST_2014_clean_resamp.csv", h=T)

heron_2014<-data.frame(Date=heron_2014$Date, Time=heron_2014$time,
                       Latitude=heron_2014$Latitude, Longitude=heron_2014$Longitude, 
                       TrackID=paste(heron_2014$Nest_ID, heron_2014$Bird_ID, sep="_"), Colony="Heron",
                       DateTime=heron_2014$DateTime, TrackTime=heron_2014$TrackTime, 
                       Vel=heron_2014$Vel, DateGMT=heron_2014$DateGMT, TimeGMT=heron_2014$TimeGMT)

results<-rbind(results, heron_2014)
write.csv(results, "GPS_2014_15_clean_resamp.csv", quote=F, row.names=F)


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

write.csv(dat_tripsplit, "GPS_2014_15_clean_resamp.csv", quote=F, row.names=F)


final_out$trip_type<-"L"
final_out[final_out$trip_id %in% trip_distances[trip_distances$time<24*4,]$trip,]$trip_type<-"S"
library(ggplot2)
qplot(data=final_out, x=Longitude, y=Latitude, colour=trip_type, geom="point")
#hmm there look to be some erroneous
strange_id<-unique(final_out[with(final_out, trip_type=="S" & Returns=="N"),]$trip_id)
strange_id2<-trip_distances[trip_distances$trip %in% strange_id,]$trip[which(trip_distances[trip_distances$trip %in% strange_id,]$fin_dist>202)]

final_out[with(final_out, trip_id %in% strange_id2),]$trip_type<-"L"
qplot(data=final_out, x=Longitude, y=Latitude, colour=trip_type, geom="point")
#better


final_out$Month<-substr(final_out$DateTime, 6,7)
final_out$Year<-substr(final_out$DateTime, 1,4)
final_out$Hr<-substr(final_out$DateTime, 12,13)
final_out$DayNight<-"Day"
final_out[final_out$Hr %in% c("19","20","21","22","23","00","01","02","03","04","05"),]$DayNight<-"Night"

write.csv(final_out, "D:/research/phd/results/GPS_2014_15_trip_fpt.csv", quote=F, row.names=F)

