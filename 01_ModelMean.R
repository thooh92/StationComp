## Author: Thomas Ohnemus
## Date: 20/02/2024
## Model Mean Calculation

# Prepare Script
rm(list = ls())

library(raster)

sub_rcp <- c("rcp45", "rcp85")
sub_var <- c("prAdjust", "tasAdjust")

# Load Data
setwd("S:/MIRO_Scenarios/interface_avail/data_weekly")
files    <- list.files(pattern = ".gri$")

# Calculate Model Mean
for(r in 1:2){ # Loop through RCPs
  for(v in 1:2){# Loop through variables
    print(paste0(sub_rcp[r], " | ", sub_var[v]))
    pat <- files[files %like% sub_rcp[r] & files %like% sub_var[v] & !files %like% "Model_Mean"]
    
    dummy <- stack(pat[1])
    s     <- stack(pat)
    
      for(i in 1:nlayers(dummy)){ 
      dummy[[i]] <- mean(s[[seq(i,nlayers(s),nlayers(dummy))]], na.rm = T) 
    # calculates model means for each dataset
    # model means assigned to dummy
    # plot(dummy[[i]], main = paste(sub_RCP[r],sub_var[v],i))
      }
    writeRaster(dummy, paste0(sub_rcp[r],"_",sub_var[v],"_Model_Mean.grd"), format = "raster", overwrite = T)
    }
}
