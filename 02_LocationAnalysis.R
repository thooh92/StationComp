## Author: Thomas Ohnemus
## Date: 30/01/2024
## Extracting and analyzing data for certain locations

# Prepare Script
rm(list = ls())

library(terra)
library(data.table)
library(raster)
library(ggplot2)

# Load data
  # Locations
locs <- data.frame(
  lat = c(51.8454),
  lon = c(10.7686)
)

locs <- vect(locs, geom=c("lon", "lat"), crs="WGS84", keepgeom=FALSE)
locs <- as(locs, "Spatial")
#plot(locs)
  
  # Climate Data
setwd("S:/MIRO_Scenarios/interface_avail/data_weekly")
files    <- list.files(pattern = ".gri$")


dates     <- seq(as.Date("1996-01-01"), as.Date("2025-12-31"), "day")
weeks     <- format(dates, format = "%G-%V")


for(i in 1:length(files)){
  print(i)
  ra <- stack(files[i])
  locs      <- spTransform(locs, CRSobj = crs(ra))
  
  ext_df    <- as.data.frame(t(raster::extract(ra, locs)))
  ext_df$w  <- unique(weeks)
  ext_df$v  <- ifelse(files[i] %like% "prAdjust", "P", "T")
  ext_df$m  <- files[i]
  
  if(i == 1){
    ext_df2 <- ext_df
  } else {
    ext_df2 <- rbind(ext_df2, ext_df)
  }
}

write.csv(ext_df2, "S:/MIRO_Scenarios/interface_avail/models_weekly.csv")

