#read in, select and compile hdf data for kernal-based extraction

rm(list=ls())
setwd("~/grive/phd/sourced_data/env_data/climatologies")
library(gdalUtils)
library(raster)
hdfz<-list.files()

sds <- get_subdatasets("TNCwwwwSmday_19951230000000_x0_X360_y-90_Y90_nx2147483647_ny2147483647.hdf")
# Any sds can then be read directly using the raster function
r <- raster(sds[1])

r1<-raster("TNCwwwwSmday_19951230000000_x0_X360_y-90_Y90_nx2147483647_ny2147483647.hdf")