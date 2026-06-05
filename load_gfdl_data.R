# Manage GFDL data sources
# Pull data from ESGF using OPeNDAP. 
# Save locally as parquet files.

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, arrow, dplyr, here) # packages (not sure if need all)

source("exposure_functions.R") # need for get_gfdl function

# make lat-lon box for the GOA
lonrange <- c(188, 235)
latrange <- c(52, 62)

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
# ph_hist_2010_2014 <- tidync(ph_hist_2010_2014_url)
# ph_hist_1990_2009 <- tidync(ph_hist_1990_2009_url)
ph_historical_files <- list(ph_hist_2010_2014_url, ph_hist_1990_2009_url, ph_ssp245_2015_3034_url)

# run functions to put GFDL pH data into a data frame for correct coords and surface layer, and compute means by month and year
ph_ssp585_surface <- get_gfdl("PH", ph_ssp585_files, "ssp585", "surface")
ph_historical_surface <- get_gfdl("PH", ph_historical_files, "historical", "surface")

# write parquet files into folder
write_parquet(ph_ssp585_surface, "data/pH/ph_ssp585_surface.parquet")
write_parquet(ph_historical_surface, "data/pH/ph_historical_surface_2020.parquet")

########## OXYGEN CONCENTRATION #########
# load future (ssp585) data for o2
o2_ssp585_2015_2034_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_201501-203412.nc"
o2_ssp585_2035_2054_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_203501-205412.nc"
o2_ssp585_2055_2074_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp585_r1i1p1f1_gr_205501-207412.nc"
# o2_ssp585_2015_2034 <- tidync(o2_ssp585_2015_2034_url)
# o2_ssp585_2035_2054 <- tidync(o2_ssp585_2035_2054_url)
# o2_ssp585_2055_2074 <- tidync(o2_ssp585_2055_2074_url)
o2_ssp585_files <- list(o2_ssp585_2015_2034_url, o2_ssp585_2035_2054_url, o2_ssp585_2055_2074_url)

# load historical data for 02
o2_hist_2010_2014_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/o2/gr/v20190726/o2_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_201001-201412.nc"
o2_hist_1990_2009_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/o2/gr/v20190726/o2_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_199001-200912.nc"
# load data from 2015-2020 from scenario ssp245
o2_ssp245_2015_2034_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Omon/o2/gr/v20180701/o2_Omon_GFDL-ESM4_ssp245_r1i1p1f1_gr_201501-203412.nc"
# o2_hist_2010_2014 <- tidync(o2_hist_2010_2014_url)
# o2_hist_1990_2009 <- tidync(o2_hist_1990_2009_url)
o2_historical_files <- list(o2_hist_2010_2014_url, o2_hist_1990_2009_url, o2_ssp245_2015_2034_url)

o2_ssp585_surface <- get_gfdl("O2", o2_ssp585_files, "ssp585", "surface")
o2_historical_surface <- get_gfdl("O2", o2_historical_files, "historical", "surface")

write_parquet(o2_ssp585_surface, "data/o2/o2_ssp585_surface.parquet")
write_parquet(o2_historical_surface, "data/o2/o2_historical_surface_2020.parquet")

########## AIR TEMPERATURE #########
# load future data
airtemp_ssp585_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Amon/tas/gr1/v20180701/tas_Amon_GFDL-ESM4_ssp585_r1i1p1f1_gr1_201501-210012.nc"
# airtemp_ssp585 <- tidync(airtemp_ssp585_url)
# load historical data
airtemp_historical_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Amon/tas/gr1/v20190726/tas_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_195001-201412.nc"
# add ssp245 for 2015-2100 to fill in gap from 2015-2020 for historical data
airtemp_ssp245_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Amon/tas/gr1/v20180701/tas_Amon_GFDL-ESM4_ssp245_r1i1p1f1_gr1_201501-210012.nc"
airtemp_historical_files <- list(airtemp_historical_url, airtemp_ssp245_url) 
# airtemp_historical <- tidync(airtemp_historical_url) # 1950-2014

airtemp_ssp585_na <- get_gfdl("AT", airtemp_ssp585_url, "ssp585", "na")
airtemp_historical_na <- get_gfdl("AT", airtemp_historical_files, "ssp585", "na")

write_parquet(airtemp_ssp585_na, "data/tas/airtemp_ssp585_na.parquet")
write_parquet(airtemp_historical_na, "data/tas/airtemp_historical_na_2020.parquet")

# PRECIPITATION
precip_ssp585_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Amon/pr/gr1/v20180701/pr_Amon_GFDL-ESM4_ssp585_r1i1p1f1_gr1_201501-210012.nc"
# precip_ssp585 <- tidync(precip_ssp585_url) # 2015-2100
precip_historical_url <- "http://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Amon/pr/gr1/v20190726/pr_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_195001-201412.nc"
# add ssp245 for 2015-2100 to fill in gap from 2015-2020 for historical data
precip_ssp245_url <- "http://esgf-node.ornl.gov/thredds/dodsC/css03_data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp245/r1i1p1f1/Amon/pr/gr1/v20180701/pr_Amon_GFDL-ESM4_ssp245_r1i1p1f1_gr1_201501-210012.nc"
precip_historical_files <- list(precip_historical_url, precip_ssp245_url)
# precip_historical <- tidync(precip_historical_url) # 1950-2014

precip_ssp585_na <- get_gfdl("PR", precip_ssp585_url, "ssp585", "na")
precip_historical_na <- get_gfdl("PR", precip_historical_files, "ssp585", "na")

write_parquet(precip_ssp585_na, "data/pr/precip_ssp585_na.parquet")
write_parquet(precip_historical_na, "data/pr/precip_historical_na_2020.parquet")


