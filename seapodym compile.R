# read in weekly seapodym nc files, compile into monthly averages and sd products
rm(list=ls())
setwd("~/grive/phd/sourced_data/SEAPODYM")
library(raster)
library(ncdf4)
library(rasterVis)
library(animation)

skj<-list.files("skj_interim_1deg_ref2015/SKJ/run-1x30d")

#use ncdf4 package to get some info
nc_open(paste("skj_interim_1deg_ref2015/SKJ/run-1x30d/", skj[1], sep=""))
# ok so lets compile adult, juvenile and total

ras<-raster(paste("skj_interim_1deg_ref2015/SKJ/run-1x30d/", skj[1], sep=""),
            varname="skj_adu")
plot(ras)

# change there 2 parameters for running on the 3 species
sp<-"yft"
sp_filepath<-"yft_interim_1deg_ref2015/YFT/run-1x30d/"

sp_filez<-list.files(sp_filepath)

monthz<-c("01", "02", "03", "04", "05", "06", "07", "08", "09",
          "10", "11", "12")

for(i in monthz)
{
  mo_filez<-sp_filez[which(substr(sp_filez,38 ,39)==i)]
  
  # adult fish potential biomass
  for(j in 1:length(mo_filez))
  {
    ras<-raster(paste(sp_filepath, mo_filez[j], sep=""),
                varname=paste(sp, "_adu", sep=""))
    if(j==1){mo_comb<-ras}else{
      mo_comb<-stack(mo_comb, ras)}
    }
  mean_mo<-calc(mo_comb, fun=mean)
  sd_mo<-calc(mo_comb, fun=sd)
  writeRaster(mean_mo, paste("month_ave_intermin_1deg_ref2015/",
                             sp, "_adu_", i, "_ave.tif", sep=""), overwrite=T)
  writeRaster(sd_mo, paste("month_ave_intermin_1deg_ref2015/",
                             sp, "_adu_", i, "_sd.tif", sep=""), overwrite=T) 
  
  # create for animation
  #if(i=="01"){adu_stack<-mean_mo}else{
   # adu_stack<-stack(adu_stack, mean_mo)}
   
  
  # juvenile fish potential biomass
  for(j in 1:length(mo_filez))
  {
    ras<-raster(paste(sp_filepath, mo_filez[j], sep=""),
                varname=paste(sp, "_juv", sep=""))
    if(j==1){mo_comb<-ras}else{
      mo_comb<-stack(mo_comb, ras)}
  }
  mean_mo<-calc(mo_comb, fun=mean)
  sd_mo<-calc(mo_comb, fun=sd)
  writeRaster(mean_mo, paste("month_ave_intermin_1deg_ref2015/",
                             sp, "_juv_", i, "_ave.tif", sep=""), overwrite=T)
  writeRaster(sd_mo, paste("month_ave_intermin_1deg_ref2015/",
                           sp, "_juv_", i, "_sd.tif", sep=""), overwrite=T)  
  
  # create for animation
  #if(i=="01"){juv_stack<-mean_mo}else{
  #  juv_stack<-stack(juv_stack, mean_mo)}
  
  
  # total fish potential biomass
  for(j in 1:length(mo_filez))
  {
    ras<-raster(paste(sp_filepath, mo_filez[j], sep=""),
                varname=paste(sp, "_tot", sep=""))
    if(j==1){mo_comb<-ras}else{
      mo_comb<-stack(mo_comb, ras)}
  }
  mean_mo<-calc(mo_comb, fun=mean)
  sd_mo<-calc(mo_comb, fun=sd)
  writeRaster(mean_mo, paste("month_ave_intermin_1deg_ref2015/",
                             sp, "_tot_", i, "_ave.tif", sep=""), overwrite=T)
  writeRaster(sd_mo, paste("month_ave_intermin_1deg_ref2015/",
                           sp, "_tot_", i, "_sd.tif", sep=""), overwrite=T)
  
  # create for animation
  #if(i=="01"){tot_stack<-mean_mo}else{
  #  tot_stack<-stack(tot_stack, mean_mo)}
  
  print(i)
}

saveGIF({
  for (i in 1:12) plot(adu_stack, i, main=i)
}, interval=1, "adu_skj.gif")


saveGIF({
  for (i in c("01", "02", "03", "04", "05",
              "06", "07", "08", "09", "10", "11", "12")){
    ras<-raster(paste("month_ave_intermin_1deg_ref2015/bet_juv_",
                      i, "_ave.tif", sep=""))
    #plot(levelplot(ras, contour=T, margin=F, main=i,at=c(0,0.01, 0.02, 0.03, 0.04,0.05, 0.06, 0.07, 0.08,0.09, 0.1, 0.15, 0.2, 0.25, 0.3)))
    plot(levelplot(ras, contour=T, margin=F, main=i))
    
  } 
}, interval=1, "bet_juv.gif")

yft<-list.files("yft_interim_1deg_ref2015/YFT/run-1x30d")

saveGIF({
  for (i in yft){
    ras<-raster(paste("yft_interim_1deg_ref2015/YFT/run-1x30d/",
                      i, sep=""), varname="yft_adu")
    #plot(levelplot(ras, contour=T, margin=F, main=i,at=c(0,0.01, 0.02, 0.03, 0.04,0.05, 0.06, 0.07, 0.08,0.09, 0.1, 0.15, 0.2, 0.25, 0.3)))
    plot(levelplot(ras, contour=T, margin=F, main=substr(i, 34, 41)))
    
  } 
}, interval=1, "yft_adu_all.gif")
#cool cool

