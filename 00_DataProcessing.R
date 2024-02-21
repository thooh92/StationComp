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
library(data.table)

setwd("S:/MIRO_Scenarios/interface_avail/data_raw") 


## Precipitation
for(g in 1:5){ # loops through GCMs
  for(q in 1:3){ # loops through RCM
    for(r in 1:2){ # loops through RCP scenarios
      for(va in 1:2){ # loops through variables of interest
        
        # nc_list loads all files based on the subsets defined earlier
        nc_list    <- intersect(list.files(pattern = sub_RCP[r]), 
                                intersect(list.files(pattern = sub_var[va]), 
                                          intersect(list.files(pattern = sub_GCM[g]), 
                                                    list.files(pattern = sub_RCM[q]))))

        print(sub_RCM[q]);print(sub_RCP[r]);print(sub_var[va]);print(length(nc_list));print("---------------")
        
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
            

            ## Produce brick from array
            b <- brick(array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), 
                       crs = CRS(st_crs(4326)$proj4string)) # using lat and lon data and correct projection to produce brick
            
            
            if(i == 1){ 
              daily  <- b 
            }
            
            
            if(i > 1){ # Here, the raster stacks of the months of respective time periods are stacked
              daily <- stack(daily, b)
            }
          }

          # Calendar Representations differ between gregorian, 365 (IPSL) and 360 (MOHC) days
          if(g == 5){
            # Extend 360 day calendar to 365 days using interpolation (following LOCA)
            v <- c()
            for(i in 1:150){  # Filling one raster randomly within a 72 day block over 30 years
              v[i] <- sample(1+(72*(i-1)):(72*i),1)
            }
            v  <- c(v, 10800)  # adding length of ra stack as last element to add rasters after last fill to stack
            
            # Calculating Mean Rasters
            # Create stack of mean rasters
            mean_rr          <- daily[[1:150]]
            
            for(i in 1:150){
              print(paste("Mean =", i))
              mean_rr[[i]]     <- mean(daily[[v[i]-1]], daily[[v[i]]], na.rm = T) 
            }
          
            # Create stack of rasters with mean in correct positions
            filled_rr <- daily[[1:v[1]]]  # Initial rasters, until first filling
            
            for(i in 1:(length(v)-1)){  # Loop filling of rasters
              print(i)
              filled_rr <- addLayer(filled_rr, mean_rr[[i]])
              filled_rr <- addLayer(filled_rr, daily[[(v[i]+1):v[i+1]]])
            }
          }
          
          if(g >= 4){  # Fill leap days for 365 day representations
            dates         <- as.Date(1:10958, origin = "1996-01-01")
            leap_days     <- which(dates %like% "-02-29")  # indices of leap days
            leap_days     <- c(leap_days, 10950)
            
            mean_leap     <- filled_rr[[1:8]]
            for(i in 1:8){
              print(paste("Mean =", i))
              mean_leap[[i]]     <- mean(filled_rr[[leap_days[i]-1]], filled_rr[[leap_days[i]]], na.rm = T) 
            }
            
            daily         <- filled_rr[[1:leap_days[1]]]
            for(i in 1:(length(leap_days)-1)){
              print(i)
              daily <- addLayer(daily, mean_leap[[i]])
              daily <- addLayer(daily, filled_rr[[(leap_days[i]+1):leap_days[i+1]]])
            }
          }

        # Summarize weekly data
        e <- list()

        
        dates <- as.Date(1:10958, origin = "1996-01-01")
        weeks <- format(dates, format = "%G-%V")
        
        for(k in 1:length(unique(weeks))){
          print(k)
          b_w   <- daily[[which(weeks == unique(weeks)[k])]]
          
          if(va == 2){
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
        if(va == 1){
          b  <- b*86400*365.00/52 # transforms precipitation unit: kg/m²*s to mm/month
        } else {
          b  <- b-273.15  # transforms temperature unit: K to °C
        }
        
        
        ## Assign correct CRS with a trick
        b <- projectRaster(b, crs = CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180 
                      +to_meter=0.01746329"))
        raster::crs(b) <- NA # Required, noone knows why though
        raster::crs(b) <- CRS(st_crs(4326)$proj4string) # sets correct Coordinate Reference System
        
        
        
                  # writes relevant data produced and assigns correct simulation names automatically
        raster::writeRaster(b, paste0("S:/MIRO_Scenarios/interface_avail/data_weekly/",
                                            sub_RCP[r],"_",sub_var[va],"_",sub_GCM[g],"_",sub_RCM[q],
                                            "_weekly_1996-2025.grd"), format = "raster", overwrite = T)
          
        gc()
        
        }}}}} 

