# 13/09/16 Lancaster, UK
# Clean and compile all GPS tracking data for use in multiple PhD chapters

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

write.csv(dat_tripsplit, "GPS_2014_15_clean_resamp_tripsplit.csv", quote=F, row.names=F)

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

#####

#data2012<-data2012[data2012$step>0,]
#data2012<-data2012[!is.na(data2012$ID),]

plot(data2012)

##experimental 4 parameter
mu0 <- c(0.0001,0.01, 0.03, 0.03)
sigma0 <- c(0.0001,0.01,0.03, 0.01)
zeromass0 <- c(0.01,0.01, 0.0001, 0.0001)
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0,pi,pi, 0)
kappa0 <- c(2,1,1, 2)
anglePar0 <- c(angleMean0,kappa0)

# parameters set up for 2 second interval data
# CURRENT runs with:2012, 2014, 2015
mu0 <- c(0.0001,0.01,0.03)
sigma0 <- c(0.0001,0.01,0.03)
zeromass0 <- c(0.01,0.0001,0.0001)
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0,pi,0)
kappa0 <- c(2,1,2)
anglePar0 <- c(angleMean0,kappa0)

m3zm <- fitHMM(data=data2sec[data2sec$year=="2012",],nbStates=3,stepPar0=stepPar0,
               anglePar0=anglePar0)

out<-data2sec[data2sec$year=="2012",]
states<-viterbi(m3zm)
out$states<-states
write.csv(out, "FINALdata2012.csv", quote=F, row.names=F)


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



## not run, need to choose HMM method then run

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVERT TRIPS OBJECT TO TRAJECTORY OBJECT TO CALCULATE TURNING ANGLES ETC.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 1. project_data
library(geosphere)

mid_point<-data.frame(centroid(cbind(final_out$Longitude, final_out$Latitude)))
final_out.Wgs <- SpatialPoints(data.frame(final_out$Longitude, final_out$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
final_out.Projected <- spTransform(final_out.Wgs, CRS=DgProj)
#now we have projected coords for WHOLE dataset (rather than cols split)

### 3. Convert to LTRAJ
library(adehabitat)

trajectories <- as.ltraj(xy=data.frame(final_out.Projected@coords[,1],
                                       final_out.Projected@coords[,2]), date=as.POSIXct(final_out$TrackTime, origin="1970/01/01", tz="GMT"), id=final_out$trip_id, typeII = TRUE)   


tt<-summary.ltraj(trajectories)
trips<-tt$id
head(trajectories[[1]])
#plot(trajectories)

## I havent run these data data cleaning loops - think they only affect the odd point

### 4. Enumerate bogus data [>3 hr time lapse]
#nonsense<-data.frame(trips, oddlocs=0)
#for (t in 1:length(trips)){
#x<-trajectories[[t]]
#nonsense$oddlocs[nonsense$trips==trips[t]]<-dim(x[x$dt>10000,])[1]
#}
#nonsense

### 5. Eliminate bogus data [>3 hr time lapse]
### removes lines with a long time interval and NA for angle

#for (t in 1:length(trips)){
#x<-trajectories[[t]]
#x<-x[!(x$dt>10000),]
#x<-x[!is.na(x$rel.angle),]
#trajectories[[t]]<-x
#}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT SINGLE HMM TO ALL TRIPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#triplist<-list()					## causes error message in viterbi: "REAL() can only be applied to a 'numeric', not a 'list'
triplist<-data.frame()
for (t in 1:length(trips)){
  uil<-trajectories[[t]]
  uil<-uil[,c(6:7,10)]
  uil$speed<-uil$dist/uil$dt
  uil$dist<-NULL
  uil$dt<-NULL
  uil[is.na(uil)]<-0
  uil$ID<-trips[[t]]
  #triplist[[t]] <- uil
  triplist <- rbind(triplist,uil)
}

hmm.2<-HMMFit(triplist[,1:2],nStates=3)		## fits a hidden markov model with 3 states - foraging and commuting and sleeping?
states<-viterbi(hmm.2,triplist[,1:2])		## extracts the predicted states for each location fix using Viterbi's algorith from the HMM fit
triplist$state2<-states$states

final_out$hmm_all_trips<-0
for(i in trips)
{
  final_out[final_out$trip_id==i,]$hmm_all_trips<-triplist[triplist$ID==i,]$state2
} ## loop to match up correct trips in triplist and final_out



########## plot separation between foraging and commuting ##########
plot(speed~rel.angle, data=triplist, type='p',col=triplist$state, pch=triplist$state)
## Nice plot :)



### SUMMARISE MAX DIST FROM COLONY AND TRIP TRAVELLING TIME FOR EACH TRIP

trip_distances<-data.frame(trip=unique(final_out$trip_id), max_dist=0, time=0, total_dist=0, fin_dist=0,
                           Returns="na",trip_type="na", Year="na", Colony="na",
                           perc_forage=0, perc_commute=0, perc_sit=0,
                           perc_forage_day=0, perc_sit_day=0)  	### create data frame for each trip

trip_distances$Returns<-as.character(trip_distances$Returns)
trip_distances$trip_type<-as.character(trip_distances$trip_type)
trip_distances$Year<-as.character(trip_distances$Year)
trip_distances$Colony<-as.character(trip_distances$Colony)


for (i in trip_distances$trip){
  trip_distances[trip_distances$trip==i,2]<-max(final_out[final_out$trip_id==i,]$ColDist)/1000
  trip_distances[trip_distances$trip==i,3]<-(max(final_out[final_out$trip_id==i,]$TrackTime)-min(final_out[final_out$trip_id==i,]$TrackTime))/3600
  trip_distances[trip_distances$trip==i,5]<-max(final_out[final_out$trip_id==i,]$ColDist)/1000
  ## Calculate distances from one point to the next and total trip distance
  x=final_out[final_out$trip_id==i,]
  x$Dist=0
  x$Dist[1]<-x$ColDist[1]/1000				### distance to first point is assumed a straight line from the nest/colony
  for (p in 2:dim(x)[1]){
    p1<-c(x$Longitude[p-1],x$Latitude[p-1])
    p2<-c(x$Longitude[p],x$Latitude[p])
    #x$Dist[p]<-pointDistance(p1,p2, lonlat=T, allpairs=FALSE)/1000			### no longer works in geosphere
    x$Dist[p]<-distMeeus(p1,p2)/1000						### great circle distance according to Meeus, converted to km
    
  }
  trip_distances[trip_distances$trip==i,4]<-sum(x$Dist)+(x$ColDist[p]/1000)	## total trip distance is the sum of all steps plus the dist from the nest of the last location - for non return trips this will be an underestimate
  trip_distances[trip_distances$trip==i,6]<-as.character(unique(final_out[final_out$trip_id==i,]$Returns))
  trip_distances[trip_distances$trip==i,7]<-as.character(unique(final_out[final_out$trip_id==i,]$trip_type))
  trip_distances[trip_distances$trip==i,8]<-as.character(unique(final_out[final_out$trip_id==i,]$Year))
  trip_distances[trip_distances$trip==i,9]<-as.character(unique(final_out[final_out$trip_id==i,]$Colony))
  trip_distances[trip_distances$trip==i,10]<-round((nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==2,])/
                                                      nrow(final_out[final_out$trip_id==i, ]))*100)
  trip_distances[trip_distances$trip==i,11]<-round((nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==3,])/
                                                      nrow(final_out[final_out$trip_id==i, ]))*100)
  trip_distances[trip_distances$trip==i,12]<-round((nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==1,])/
                                                      nrow(final_out[final_out$trip_id==i, ]))*100)
  
  
  trip_distances[trip_distances$trip==i,13]<-round((nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==2 & final_out$DayNight=="Day",])/
                                                      nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==2,]))*100)
  trip_distances[trip_distances$trip==i,14]<-round((nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==1 & final_out$DayNight=="Day",])/
                                                      nrow(final_out[final_out$trip_id==i &final_out$hmm_all_trips==1,]))*100)
  
}

trip_distances$days<-ceiling(trip_distances$time/24)

## using merge :) to attrib days to final_out

#m1<-merge(final_out, trip_distances, by.x="trip_id",by.y="trip", set.all.x=TRUE)
#final_out$days<-m1$days ## doesnt work.. need to sort rows properly

#bodge but works
final_out$days<-0
for(i in trip_distances$trip)
{
  final_out[final_out$trip_id==i,]$days<-trip_distances[trip_distances$trip==i,]$days
}




write.csv(trip_distances, "D:/research/phd/results/GPS_2014_15_trip_summary.csv", quote=F, row.names=F)

trip_distances<-read.csv("D:/research/phd/results/GPS_2014_15_trip_summary.csv", h=T)

