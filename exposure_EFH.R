# ============================================================================
# Exposure script
# ============================================================================
# Purpose: This script calculates exposure for all species and all exposure factors in the CVA.

# clear work space and free up memory
rm(list = ls())
gc()

# load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
  pacman::p_load(ggplot2, stringr, tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, terra, readxl, data.table, gstat, ggpubr, cowplot, patchwork, heatwaveR)
    pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
    
# identify wd using here package
here::i_am("exposure_EFH.R")

# load custom analysis functions
source(here("exposure_functions.R")) 

# ============================================================================
# SECTION 1A: Load and Filter ROMS Data - Historical and Future
# ============================================================================
# The Arrow package allows "lazy evaluation" - you can filter the data before
# loading it into RAM, which is much faster than loading everything first. Load
# marine heatwave index.
# ============================================================================

# load_roms function comes from exposure_functions.R
load_roms("surface", "temp") # sea surface temp
load_roms("bottom", "temp") # bottom temp
load_roms("surface", "salt") # surface salinity 
load_roms("bottom", "salt") # bottom salinity 
load_roms("surface", "PhL") # large phytoplankton
load_roms("surface", "PhS") # small phytoplankton
load_roms("surface", "Cop") # small copepod concentration
load_roms("surface", "NCa") # large copepod concentration
load_roms("surface", "Eup") # euphausiid concentration
load_roms("surface", "MZL") # large microzooplankton concentration
load_roms("surface", "MZS") # small microzooplankton concentration

# combine phytoplankton (small and large)
PhT_ssp585_surface <- rbind(PhL_ssp585_surface, PhS_ssp585_surface) 
PhT_ssp585_surface <- PhT_ssp585_surface |> summarize(value_dc = sum(value_dc), .by = c(cell_id, lon_rho, lat_rho, year, month))
PhT_hindcast_surface <- rbind(PhL_hindcast_surface, PhS_hindcast_surface) 
PhT_hindcast_surface <- PhT_hindcast_surface |> summarize(value = sum(value), .by = c(cell_id, lon_rho, lat_rho, year, month))

# combine zooplankton (all)
ZP_ssp585_surface <- rbind(Cop_ssp585_surface, NCa_ssp585_surface, Eup_ssp585_surface, MZL_ssp585_surface, MZS_ssp585_surface) 
ZP_ssp585_surface <- ZP_ssp585_surface |> summarize(value_dc = sum(value_dc), .by = c(cell_id, lon_rho, lat_rho, year, month))
ZP_hindcast_surface <- rbind(Cop_hindcast_surface, NCa_hindcast_surface, Eup_hindcast_surface, MZL_hindcast_surface, MZS_hindcast_surface) 
ZP_hindcast_surface <- ZP_hindcast_surface |> summarize(value = sum(value), .by = c(cell_id, lon_rho, lat_rho, year, month))

# load and process mhw index (current metric: total intensity)
mhw_ssp585 <- readRDS(here("data/mhw/processed_mhw_ssp585.RDS"))
mhw_hindcast <- readRDS(here("data/mhw/processed_mhw_hindcast.RDS"))

# ============================================================================
# SECTION 1B: Load GFDL Data - Historical and Future
# ============================================================================
# Load parquet files from data folder. Plot variogram. Krige.
# These files are pulled and processed using the load_gfdl_data.R script.
# ============================================================================

# load prediction grid (10 km^2 resolution)
load(here::here("data", "10km_grid.rda"))
prediction.grid = as.data.frame(grid)
grid_sf = st_as_sf(prediction.grid, coords = c("lon_rho", "lat_rho"), crs = 4326) |>
  st_shift_longitude()

# fit and print variogram model results
# identify psill, range, and nugget
find_vgm_values <- function(vgm_data){
  
  this_variable.vgm = variogram(mean_annual ~ 1, data = vgm_data)
  plot(this_variable.vgm)
  this_variable.vgm_fit = fit.variogram(this_variable.vgm, model = vgm(psill = 0.1, model = "Gau", range = 450, nugget = 0.1))
  plot(this_variable.vgm, this_variable.vgm_fit)
  
  # print and return table with values 
  print(this_variable.vgm_fit)
  return(this_variable.vgm_fit)
}

# krige by year
# identify the best fit psill, nugget, and range values for each year 
run_krige <- function(data){
  for_kriging <- data |>
    summarize(mean_annual = mean(value), .by = c(cell_id, lon, lat, year)) |> # need only one point per location to krige by year 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  kriged_list = list()
  for(i in unique(for_kriging$year)) {
    df = subset(for_kriging, year == i)
    
    # identify year-specific psill, nugget, and range values
    this_variable.vgm_fit <- find_vgm_values(df)
    
    nugget_val <- this_variable.vgm_fit[1, "psill"] # nugget
    sill_val   <- this_variable.vgm_fit[2, "psill"] # partial sill
    range_val  <- this_variable.vgm_fit[2, "range"] # range
    
    vgm = variogram(mean_annual ~ 1, data = df)
    vgm_fit = fit.variogram(vgm, model = vgm(psill = sill_val, model = "Gau", range = range_val, nugget = nugget_val))
    
    krig = gstat::krige(mean_annual ~ 1, df, grid_sf, model = vgm_fit)
    
    krig = as.data.frame(krig)
    krig$year = i
    kriged_list[[i]] = krig 
  }
  all_kriged = dplyr::bind_rows(kriged_list) |>
    rename(kriged = var1.pred)
  
  return(all_kriged)
}

# load all GFDL data files and filter to correct time periods
pH_ssp585_surface <- open_dataset(here("data/pH/ph_ssp585_surface.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

pH_historical_surface <- open_dataset(here("data/pH/ph_historical_surface.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

pH_ssp585_bottom <- open_dataset(here("data/pH/ph_ssp585_bottom.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

pH_historical_bottom <- open_dataset(here("data/pH/ph_historical_bottom.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

o2_ssp585_surface <- open_dataset(here("data/o2/o2_ssp585_surface.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

o2_historical_surface <- open_dataset(here("data/o2/o2_historical_surface.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

o2_ssp585_bottom <- open_dataset(here("data/o2/o2_ssp585_bottom.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

o2_historical_bottom <- open_dataset(here("data/o2/o2_historical_bottom.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

AT_ssp585 <- open_dataset(here("data/tas/airtemp_ssp585_na.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

AT_historical <- open_dataset(here("data/tas/airtemp_historical_na_2020.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

PR_ssp585 <- open_dataset(here("data/pr/precip_ssp585_na.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
  filter(year != 2035 & year != 2054) |> # remove years with kriging issues
    collect()

PR_historical <- open_dataset(here("data/pr/precip_historical_na_2020.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

gfdl_files <- c("pH_ssp585_surface", "pH_historical_surface", "pH_ssp585_bottom", "pH_historical_bottom", "o2_ssp585_surface", "o2_historical_surface", "o2_ssp585_bottom", "o2_historical_bottom", "AT_ssp585", "AT_historical", "PR_ssp585", "PR_historical")
  
kriged <- c()
for(i in 1:length(gfdl_files)){
  this_data = get(gfdl_files[i]) 
  kriged[[gfdl_files[i]]] <- run_krige(this_data)
}

exposure_factors <- c("SST", "BT", "SS", "BS", "PhT", "ZP", "SpH", "BpH", "SO2", "BO2", "AT", "PR", "MHW") 
# sea-surface temp, bottom temp, surface salinity, bottom salinity,
# phytoplankton concentration, zooplankton concentration, pH (ocean acidification),
# oxygen concentration, air temperature, precipitation, mhw index

# ============================================================================
# SECTION 2A: Calculate Anomalies and Variabilities
# ============================================================================
# Calculate standardized anomalies using a Z-score: (future mean - historical mean) / historical SD
# for each exposure factor. Calculate variability using an F-test: (future variance / 
# historical variance) for each exposure factor. 
# ============================================================================

SST <- create_anomaly_roms(temp_ssp585_surface, temp_hindcast_surface) |>
  mutate(exposure_factor = "SST")

BT <- create_anomaly_roms(temp_ssp585_bottom, temp_hindcast_bottom) |>
  mutate(exposure_factor = "BT")

SS <- create_anomaly_roms(salt_ssp585_surface, salt_hindcast_surface) |>
  mutate(exposure_factor = "SS")

BS <- create_anomaly_roms(salt_ssp585_bottom, salt_hindcast_bottom) |>
  mutate(exposure_factor = "BS")

PhT <- create_anomaly_roms(PhT_ssp585_surface, PhT_hindcast_surface) |>
    mutate(exposure_factor = "PhT")

ZP <- create_anomaly_roms(ZP_ssp585_surface, ZP_hindcast_surface) |>
    mutate(exposure_factor = "ZP")

# calculate anomaly for ESM variables (note different function in exposure_functions.R)

SpH <- create_anomaly_gfdl("SpH", kriged[["pH_ssp585_surface"]], kriged[["pH_historical_surface"]]) |>
    mutate(exposure_factor = "SpH")

BpH <- create_anomaly_gfdl("BpH", kriged[["pH_ssp585_bottom"]], kriged[["pH_historical_bottom"]]) |>
    mutate(exposure_factor = "BpH")

SO2 <- create_anomaly_gfdl("SO2", kriged[["o2_ssp585_surface"]], kriged[["o2_historical_surface"]]) |>
    mutate(exposure_factor = "SO2")

BO2 <- create_anomaly_gfdl("BO2", kriged[["o2_ssp585_bottom"]], kriged[["o2_historical_bottom"]]) |>
    mutate(exposure_factor = "BO2")

AT <- create_anomaly_gfdl("AT", kriged[["AT_ssp585"]], kriged[["AT_historical"]]) |>
    mutate(exposure_factor = "AT")

PR <- create_anomaly_gfdl("PR", kriged[["PR_ssp585"]], kriged[["PR_historical"]]) |>
    mutate(exposure_factor = "PR")
  
# calculate anomaly for chosen MHW index 
MHW <- create_anomaly_mhw(mhw_ssp585, mhw_hindcast) |>
  mutate(exposure_factor = "MHW")

# list of exposure factor data frames
exposure_factors_df <- list(SST, BT, SS, BS, PhT, ZP, SpH, BpH, SO2, BO2, AT, PR, MHW)
names(exposure_factors_df) <- c("SST", "BT", "SS", "BS", "PhT", "ZP", "SpH", "BpH", "SO2", "BO2", "AT", "PR", "MHW")

sf::sf_use_s2(FALSE) # IMPORTANT - turns off spherical geometry

# for plotting:
# load coastline 
coast <- ne_coastline(scale = "medium", returnclass = "sf") |>
  st_crop(xmin = -172, xmax = -130, ymin = 50, ymax = 62) |>  # crop to GOA region
    st_shift_longitude()  # convert to 0-360° longitude to match ROMS data

# load GOA shading
GOA = ne_countries(scale = "medium", 
                   returnclass = "sf") |>
  st_crop(xmin = -172, xmax = -130, ymin = 50, ymax = 62) |>  
    st_shift_longitude()  

# load Canada shading
canada = ne_countries(scale = "medium", country = "Canada", returnclass = "sf") |>
  st_crop(xmin = -172, xmax = -130, ymin = 50, ymax = 62) |>  
    st_shift_longitude()  

# create single anomaly plots for all exposure factors and save locally (13 plots)
for (k in 1:length(exposure_factors_df)){
  create_anomaly_plot(exposure_factors_df[[k]], exposure_factors[k])
}

# create series of plots: show average mean, sd, future change, and anomaly (13 plots)
for (k in 1:length(exposure_factors_df)){
  create_all_anomaly_plots(exposure_factors_df[[k]], exposure_factors[k])
}

# create series of plots for variability: future variance, historical variance, and variability (13 plots)
for (k in 1:length(exposure_factors_df)){
  create_variability_plots(exposure_factors_df[[k]], exposure_factors[k])
}

# ============================================================================
# SECTION 3: Load and plot EFH maps.
# ============================================================================
# Source and create EFH maps for each species by transforming .gdb files into 
# sf objects. 
# ============================================================================

# file paths for species distribution data 
gdb_path <- here::here("data", "GOA_groundfish_2023.gdb")
scallop_path <- here::here("data", "EFH_2018_Scallop.gdb")
bts_sdm_path <- here::here("data", "bts_sdms/")
diet_derived_path <- here::here("data", "diet_sdms/")
depth_temp_path <- here::here("data", "depth_temp_sdms/")

# list of all species, associated distribution file/layer, and exposure assignments
layer_names <- read_excel(here("exposure_assignments.xlsx"), sheet = "CVA_exposure_assignments", skip = 1,
                           col_names = c("group", "species_name", "path", "EFH_level", "layer", "sea_temperature", 
                                         "salinity", "ph", "phytoplankton", "zooplankton",
                                         "oxygen", "air_temperature", "precipitation", "marine_heatwave", "data_quality")) |>
  drop_na(layer)

# clean layer_names and make lists
layer_names$air_temperature[layer_names$air_temperature == "NA"] <- NA
layer_names$precipitation[layer_names$precipitation == "NA"] <- NA
species_layers <- layer_names$layer
species_name <- layer_names$species_name
paths <- layer_names$path
EFH_level <- layer_names$EFH_level

# # test distribution plots
# plot_species_distribution("Redbanded rockfish")
# plot_species_distribution("Weathervane scallop")
# plot_species_distribution("Red king crab")
# plot_species_distribution(Spot shrimp")
# plot_species_distribution("Geoduck clam")

# for loop to create EFH plots for all species (57 plots)
for (i in 1:length(species_layers)) {
  this_species_name <- species_name[i]
  plot_species_distribution(this_species_name)
}

# ============================================================================
# SECTION 4: Calculate exposure scores.
# ============================================================================
# Create overlap of EFH and ROMS outputs. Calculate exposure scores for anomaly
# and variability exposure factors.
# ============================================================================

# test overlap function and calculating exposure scores
geoduck_clam <- create_overlap("Geoduck clam", exposure_factors_df[["MHW"]])
calculate_exposure_score(geoduck_clam, "anomaly")

# create data frame to place exposure scores in 
exposure_scores <- layer_names |>
  dplyr::select(!c(path, EFH_level, layer)) |>
  pivot_longer(cols=c("sea_temperature", "salinity", "ph", "phytoplankton", "zooplankton",
                      "oxygen", "air_temperature", "precipitation", "marine_heatwave"),
                               names_to = "variable",
                               values_to = "exposure_factor") |>
  dplyr::select(!variable) |>
    drop_na(exposure_factor) 

# calculate exposure scores for all species and anomalies and place them in the df
m <- 1
anomaly_scores <- exposure_scores
variability_scores <- exposure_scores
for (i in 1:length(species_layers)) {
  # identify appropriate list of exposure factors for each species
  exposure_factors_list <- layer_names[i,6:14]
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])

  for(k in 1:length(these_exposure_factors$V1)){
    original_exposure_data <- create_overlap(species_name[i], exposure_factors_df[[paste(these_exposure_factors$V1[k])]], these_exposure_factors$V1[k])
    anomaly_score <- calculate_exposure_score(original_exposure_data, "anomaly")
    variability_score <- calculate_exposure_score(original_exposure_data, "variability")
    
    anomaly_scores[m,"score"] <- anomaly_score
    anomaly_scores[m,"type"] <- "anomaly"
    
    variability_scores[m,"score"] <- variability_score
    variability_scores[m,"type"] <- "variability"

    m <- m + 1
  }
}

# remove variability for MHW index
all_exposure_scores <- rbind(anomaly_scores, variability_scores) |>
  filter(!(exposure_factor == "MHW" & type == "variability"))

# save csv with exposure scores
write.csv(all_exposure_scores, "results/exposure_factor_scores_072926.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: Summary plots
# ============================================================================

# ============================================================================

# plot distribution of scores for all exposure factors for each species
  # test plots
  exposure_histogram_series("Walleye pollock", "anomaly")
  exposure_histogram_series("Walleye pollock", "variability")
  
  # create series of plots for all species
  # separate files for anomaly and variability
  # for (i in 1:length(species_layers)) {
  #   exposure_histogram_series(species_name[i], "anomaly")
  #   exposure_histogram_series(species_name[i], "variability")
  # }

# # plot:
  # # 1. exposure map (binned anomalies within species distribution)
  # # 2. histogram of anomalies across the species distribution
  # # 3. distribution of exposure scores and weighted average 
  
  # # for loop creates 3 exposure plots for all EFH species and assigned exposure factors (warning: creates ~900 plots)
  # for (i in 1:length(species_layers)) {
  #   # identify appropriate list of exposure factors for each species
  #   exposure_factors_list <- layer_names[i,6:13]
  #   these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  #   these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])
  #   
  #   for(k in 1:length(these_exposure_factors$V1)){
  #     create_exposure_plots(species_name[i], exposure_factors_df[[paste(these_exposure_factors$V1[k])]], these_exposure_factors$V1[k])
  #   }
  # }

# plot: exposure map with inset anomaly/variability histogram
  # saves 12-17 plots for each species to a folder -- to combine manually

# focal species plots
group_exposure_plots("Walleye pollock", "anomaly")
group_exposure_plots("Walleye pollock", "variability")
group_exposure_plots("Pacific cod", "anomaly")
group_exposure_plots("Pacific cod", "variability")
group_exposure_plots("Petrale sole", "anomaly")
group_exposure_plots("Petrale sole", "variability")
group_exposure_plots("Shortraker rockfish", "anomaly")
group_exposure_plots("Shortraker rockfish", "variability")
group_exposure_plots("Razor clam", "anomaly")
group_exposure_plots("Razor clam", "variability")
group_exposure_plots("Pacific herring", "anomaly")
group_exposure_plots("Pacific herring", "variability")
group_exposure_plots("Sockeye salmon", "anomaly")
group_exposure_plots("Sockeye salmon", "variability")
group_exposure_plots("Pacific spiny dogfish", "anomaly")
group_exposure_plots("Pacific spiny dogfish", "variability")