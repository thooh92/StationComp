## Author: Thomas Ohnemus
## Date: 30/01/2024
## Processing Climate Projection Data to find best Model for Location

# Prepare Script
write("TMPDIR = 'S:/temp'", file =file.path(Sys.getenv('R_USER'), '.Renviron'))
rm(list = ls())

# Subsets
sub_GCM <- c("MPI-M-MPI-ESM-LR", "CNRM-CERFACS-CNRM-CM5", "ICHEC-EC-EARTH", "IPSL-IPSL-CM5A-MR",
             "MOHC-HadGEM2-ES") # used to loop through GCMs
sub_RCM <- c("RCA4", "CCLM4-8-17", "RACMO22E") # used to loop through RCMs
sub_RCP <- c("rcp45", "rcp85") # used to loop through RCP scenarios
sub_var <- c("prAdjust", "tasAdjust")

# Packages
library(ncdf4)
library(raster)
library(sf)

setwd("S:/MIRO_Scenarios/interface_avail/data_raw") 

## Precipitation
for(g in 2:2){ # loops through GCMs
  for(q in 2:2){ # loops through RCM
    for(r in 1:2){ # loops through RCP scenarios
      for(v in 1:2){ # loops through variables of interest
        
        # nc_list loads all files based on the subsets defined earlier
        nc_list    <- intersect(list.files(pattern = sub_RCP[r]), 
                                intersect(list.files(pattern = sub_var[v]), 
                                          intersect(list.files(pattern = sub_GCM[g]), 
                                                    list.files(pattern = sub_RCM[q]))))

        print(sub_RCM[q]);print(sub_RCP[r]);print(sub_var[v]);print(length(nc_list));print("---------------")
        
        if(length(nc_list) > 0){ # if combination of GCM and RCM is nonexistent, loop jumps to next combination
          for(i in 1:length(nc_list)){ # this loops through all .nc files of a single simulation run
            
            ## Get relevant Metadata
            nc_data <- nc_open(nc_list[i]) # Grabs ith decade
            var <- sub("\\_.*", "", nc_list[i]) # Creates variable name from .nc metadata
            lon <- ncvar_get(nc_data, "rlon") # Creates lon data from .nc metadata
            lat <- ncvar_get(nc_data, "rlat", verbose = F) # Creates lat data from .nc metadata

            
            ## Produce Array from data
            array <- ncvar_get(nc_data, var) # gets data in an array
            fillvalue <- ncatt_get(nc_data, var, "_FillValue") # Extracts fill value from .nc metadata
            array[array == fillvalue$value] <- NA # fill values are replaced with NA
              
            
            #################
            year <- as.POSIXct(substr(nc_list[i], start = nchar(nc_list[i])-19,
                           stop = nchar(nc_list[1])-12), format = "%Y%m%d")
            dates <- as.Date(1:dim(array)[3], origin = year)
            weeks <- format(dates, format = "%G-%V")
            #################
 
            
            ## Produce brick from array
            b <- brick(array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), 
                       crs = CRS(st_crs(4326)$proj4string)) # using lat and lon data and correct projection to produce brick
            e <- list()
            
            # Summarize weekly data
            for(k in 1:length(unique(weeks))){
              print(k)
              b_w   <- b[[which(weeks == unique(weeks)[k])]]
              
              if(v == 2){
                b_w   <- mean(b_w, na.rm = T)
              } else {
                b_w   <- sum(b_w, na.rm = T)
              }
              
              e[[k]]<- b_w
            }
            b <- stack(e)
            names(b) <- unique(weeks)
            
            
            b <- t(b) # transpose brick
            b <- raster::flip(b, direction = "y") # flip brick, due to weird .nc data storage

            # Correcting values
            if(v == 1){
              array  <- array*86400*365.00/52 # transforms precipitation unit: kg/m²*s to mm/month
            } else {
              array  <- array-273.15  # transforms temperature unit: K to °C
            }
   
            
            
            
            ## Assign correct CRS with a trick
            b <- projectRaster(b, crs = CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180 
                      +to_meter=0.01746329"))
            raster::crs(b) <- NA # Required, noone knows why though
            raster::crs(b) <- CRS(st_crs(4326)$proj4string) # sets correct Coordinate Reference System

            
            if(i == 1){ # Creates "dummy" data in the first loop, !!AT THE MOMENT FOR ONE DECADE!!
              daily  <- b # dummy is 1800 rasters long, sum of daily data per decade
            }
            
            
            if(i > 1){ # Here, the raster stacks of the months of respective time periods are stacked
              daily <- stack(daily, b)
            }
          }
          
          # writes relevant data produced and assigns correct simulation names automatically
          raster::writeRaster(daily, paste0("S:/MIRO_Scenarios/interface_avail/data_weekly/",
                                            sub_RCP[r],"_",sub_var[v],"_",sub_GCM[g],"_",sub_RCM[q],
                                            "_weekly_1996-2025.grd"), format = "raster", overwrite = T)
          
          gc()
        }}}}} 

