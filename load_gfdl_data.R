# Manage GFDL data sources
# Pull data from ESGF using OPeNDAP. 
# Save locally as parquet files.

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, arrow, dplyr, here, data.table) # packages 
here::i_am("load_gfdl_data.R")

# make lat-lon box for the GOA
lonrange <- c(188, 235)
latrange <- c(52, 62)

# For a particular exposure factor, process GFDL data, extract necessary parts, 
# and save as a data frame for exposure analysis.
# 
# exposure_factor: Character exposure factor name.
# data: List of OPeNDAP urls (that link to NetCDF files) from historical or future time periods.
# run: Either "historical" or "ssp585" (future).
# esm_slice: "surface" or "bottom" 
get_gfdl <- function(exposure_factor, data, run, esm_slice){
  # loop over esm files
  this_esm_scenario <- list()
  for(i in 1:length(data)){
    esm_file <- data[[i]] # input list of urls to pull ESM data from
    esm_data <- tidync(esm_file)  
    
    origin <- ncmeta::nc_atts(esm_file, "time") %>% 
      mutate(across(variable:value, as.character)) %>% # fixes classes to pull origin
      tidyr::unnest(cols = c(value)) %>%
      filter(name == 'units') %>%
      dplyr::select(value) %>%
      mutate(value = str_replace(value, 'days since ', '')) %>%
      pull(value) %>%
      as.Date()
    
    # pull data from particular latitude and longitude
    this_esm_data <- esm_data %>% hyper_filter(lon = lon > lonrange[1] & lon <= lonrange[2], 
                                               lat = lat > latrange[1] & lat <= latrange[2]) 
    
    # get surface or bottom slice
    if(esm_slice == "surface"){
      this_esm_data <- this_esm_data %>%
        hyper_tibble %>% # added to make group_by work 
        group_by(time, lon, lat) %>%
        slice_min(lev) %>% # take surface slice
        ungroup()
      
      this_esm_data$time <- as.Date(this_esm_data$time)
      
      # add time and cell index
      this_esm_data <- this_esm_data %>% 
        mutate(year = year(time),
               month = month(time),
               cell_id = paste(lon, lat, sep = '_'))
    }
    else if(esm_slice == "bottom"){
      this_esm_data <- this_esm_data %>%
        hyper_tibble %>% # added to make group_by work 
        group_by(time, lon, lat) %>%
        slice_max(lev) %>% # take bottom slice
        ungroup()
      
      this_esm_data$time <- as.Date(this_esm_data$time)
      
      # add time and cell index
      this_esm_data <- this_esm_data %>% 
        mutate(year = year(time),
               month = month(time),
               cell_id = paste(lon, lat, sep = '_'))
    }
    else{ # for non-ocean variables (no surface or bottom slice needed)
      this_esm_data$time <- as.Date(this_esm_data$time)
      
      this_esm_data <- this_esm_data %>% 
        hyper_tibble %>%
        mutate(year = year(time),
               month = month(time),
               cell_id = paste(lon, lat, sep = '_'))
    }
    
    # standardize exposure factor column name
    this_esm_data <- this_esm_data %>%
      rename(any_of(c(value = "ph", value = "o2", value = "tas", value = "pr"))) |>
      mutate(run = run,
             slice = esm_slice)
    
    this_esm_scenario[[i]] <- this_esm_data
  }
  this_esm_scenario_df <- rbindlist(this_esm_scenario)
  
  return(this_esm_scenario_df)
}

########## PH (OCEAN ACIDIFICATION) #########
# load future data (scenario ssp585) for ph
ph_ssp585_2015_2034_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/ph/gr/v20180701/ph_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-203412.nc"
ph_ssp585_2035_2054_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/ph/gr/v20180701/ph_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_203501-205412.nc"
ph_ssp585_2055_2074_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/ph/gr/v20180701/ph_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_205501-207412.nc"
# ph_ssp585_2015_2034 <- tidync(ph_ssp585_2015_2034_url)
# ph_ssp585_2035_2054 <- tidync(ph_ssp585_2035_2054_url)
# ph_ssp585_2055_2074 <- tidync(ph_ssp585_2055_2074_url)
ph_ssp585_files <- list(ph_ssp585_2015_2034_url, ph_ssp585_2035_2054_url, ph_ssp585_2055_2074_url)
# ph_ssp585_files <- list(ph_ssp585_2015_2034, ph_ssp585_2035_2054, ph_ssp585_2055_2074)

# load historical data for ph
ph_hist_2010_2014_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/ph/gr/v20190726/ph_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc"
ph_hist_1990_2009_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/ph/gr/v20190726/ph_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc"

# load data from 2015-2020 from scenario ssp245
ph_ssp245_2015_3034_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Omon/ph/gr/v20180701/ph_Omon_GFDL-ESM4_ssp245_r1i1p1f1_gr_201501-203412.nc"
ph_historical_files <- list(ph_hist_2010_2014_url, ph_hist_1990_2009_url, ph_ssp245_2015_3034_url)

# run functions to put GFDL pH data into a data frame for correct coords and surface layer, and compute means by month and year
ph_ssp585_surface <- get_gfdl("PH", ph_ssp585_files, "ssp585", "surface")
ph_ssp585_bottom <- get_gfdl("PH", ph_ssp585_files, "ssp585", "bottom")
ph_historical_surface <- get_gfdl("PH", ph_historical_files, "historical", "surface")
ph_historical_bottom <- get_gfdl("PH", ph_historical_files, "historical", "bottom")

# write parquet files into folder
write_parquet(ph_ssp585_surface, "data/pH/ph_ssp585_surface.parquet")
write_parquet(ph_historical_surface, "data/pH/ph_historical_surface.parquet")
write_parquet(ph_ssp585_bottom, "data/pH/ph_ssp585_bottom.parquet")
write_parquet(ph_historical_bottom, "data/pH/ph_historical_bottom.parquet")

########## OXYGEN CONCENTRATION #########
# load future (ssp585) data for o2
o2_ssp585_2015_2034_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-203412.nc"
o2_ssp585_2035_2054_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_203501-205412.nc"
o2_ssp585_2055_2074_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_205501-207412.nc"
o2_ssp585_files <- list(o2_ssp585_2015_2034_url, o2_ssp585_2035_2054_url, o2_ssp585_2055_2074_url)

# load historical data for 02
o2_hist_2010_2014_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/o2/gr/v20190726/o2_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc"
o2_hist_1990_2009_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/o2/gr/v20190726/o2_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc"

# load data from 2015-2020 from scenario ssp245
o2_ssp245_2015_2034_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp245_r1i1p1f1_gr_201501-203412.nc"
o2_historical_files <- list(o2_hist_2010_2014_url, o2_hist_1990_2009_url, o2_ssp245_2015_2034_url)

o2_ssp585_surface <- get_gfdl("O2", o2_ssp585_files, "ssp585", "surface")
o2_historical_surface <- get_gfdl("O2", o2_historical_files, "historical", "surface")
o2_ssp585_bottom <- get_gfdl("O2", o2_ssp585_files, "ssp585", "bottom")
o2_historical_bottom <- get_gfdl("O2", o2_historical_files, "historical", "bottom")

write_parquet(o2_ssp585_surface, "data/o2/o2_ssp585_surface.parquet")
write_parquet(o2_historical_surface, "data/o2/o2_historical_surface.parquet")
write_parquet(o2_ssp585_bottom, "data/o2/o2_ssp585_bottom.parquet")
write_parquet(o2_historical_bottom, "data/o2/o2_historical_bottom.parquet")

########## AIR TEMPERATURE #########
# load future data
airtemp_ssp585_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Amon/tas/gr1/v20180701/tas_Amon_GFDL-ESM4_ssp585_r1i1p1f1_gr1_201501-210012.nc"

# load historical data
airtemp_historical_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Amon/tas/gr1/v20190726/tas_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_195001-201412.nc"

# add ssp245 for 2015-2100 to fill in gap from 2015-2020 for historical data
airtemp_ssp245_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Amon/tas/gr1/v20180701/tas_Amon_GFDL-ESM4_ssp245_r1i1p1f1_gr1_201501-210012.nc"
airtemp_historical_files <- list(airtemp_historical_url, airtemp_ssp245_url) 

airtemp_ssp585_na <- get_gfdl("AT", airtemp_ssp585_url, "ssp585", "na")
airtemp_historical_na <- get_gfdl("AT", airtemp_historical_files, "ssp585", "na")

write_parquet(airtemp_ssp585_na, "data/tas/airtemp_ssp585_na.parquet")
write_parquet(airtemp_historical_na, "data/tas/airtemp_historical_na_2020.parquet")

######### PRECIPITATION ############
precip_ssp585_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Amon/pr/gr1/v20180701/pr_Amon_GFDL-ESM4_ssp585_r1i1p1f1_gr1_201501-210012.nc"
precip_historical_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Amon/pr/gr1/v20190726/pr_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_195001-201412.nc"

# add ssp245 for 2015-2100 to fill in gap from 2015-2020 for historical data
precip_ssp245_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Amon/pr/gr1/v20180701/pr_Amon_GFDL-ESM4_ssp245_r1i1p1f1_gr1_201501-210012.nc"
precip_historical_files <- list(precip_historical_url, precip_ssp245_url)

precip_ssp585_na <- get_gfdl("PR", precip_ssp585_url, "ssp585", "na")
precip_historical_na <- get_gfdl("PR", precip_historical_files, "ssp585", "na")

write_parquet(precip_ssp585_na, "data/pr/precip_ssp585_na.parquet")
write_parquet(precip_historical_na, "data/pr/precip_historical_na_2020.parquet")


