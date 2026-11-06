# ============================================================================
# GOA CVA exposure analysis 
# ============================================================================
# Purpose: This script calculates exposure scores for all species and all 
# exposure factors in the CVA.
# Author: Isabelle Galko
# Date: November 6, 2026

# load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
  pacman::p_load(ggplot2, stringr, tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, terra, readxl, data.table, gstat, ggpubr, cowplot, patchwork, heatwaveR)
    pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
    
# identify wd using here package
here::i_am("exposure_EFH.R")

# load custom analysis functions
source(here("exposure_functions.R")) 

sf::sf_use_s2(FALSE) # necessary for sf plots - turns off spherical geometry

# ============================================================================
# SECTION 1: Load and filter ROMS data - historical and future
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
# SECTION 2: Load GFDL model outputs - historical and future 
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

# write.csv(pH_ssp585_surface, file = "SpH_ssp585.csv", row.names = FALSE)

pH_historical_surface <- open_dataset(here("data/pH/ph_historical_surface.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

# write.csv(pH_historical_surface, file = "SpH_historical.csv", row.names = FALSE)

pH_ssp585_bottom <- open_dataset(here("data/pH/ph_ssp585_bottom.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

# write.csv(pH_ssp585_bottom, file = "BpH_ssp585.csv", row.names = FALSE)

pH_historical_bottom <- open_dataset(here("data/pH/ph_historical_bottom.parquet")) |>
  filter(year >= 1993 & year <= 2020) |>
    collect()

# write.csv(pH_historical_bottom, file = "BpH_historical.csv", row.names = FALSE)

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
  filter(year >= 1993 & year <= 2020) |> # remove problem years
    collect()

gfdl_files <- c("pH_ssp585_surface", "pH_historical_surface", "pH_ssp585_bottom", "pH_historical_bottom", "o2_ssp585_surface", 
                "o2_historical_surface", "o2_ssp585_bottom", "o2_historical_bottom", "AT_ssp585", "AT_historical", "PR_ssp585", "PR_historical")
  
kriged <- c()
for(i in 1:length(gfdl_files)){
  this_data = get(gfdl_files[i]) 
  kriged[[gfdl_files[i]]] <- run_krige(this_data)
}

exposure_factors <- c("SST", "BT", "SS", "BS", "PhT", "ZP", "SpH", "BpH", "SO2", "BO2", "AT", "PR", "MHW") 
factor_full_names <- c("Sea surface temperature", "Bottom temperature", "Surface salinity", "Bottom salinity", 
                       "Phytoplankton concentration", "Zooplankton concentration",
                       "Surface pH", "Bottom pH", "Surface oxygen concentration", "Bottom oxygen concentration", 
                       "Air temperature", "PR" = "Precipitation", "MHW" = "Marine heatwaves")

# ============================================================================
# SECTION 3: Calculate anomalies and variabilities for all exposure factors
# ============================================================================
# Calculate standardized anomalies using a Z-score: (future mean - historical mean) / historical SD
# for each exposure factor. Calculate variability using an F-test: (future variance / 
# historical variance) for each exposure factor. 
# ============================================================================

SST <- create_anomaly_roms(temp_ssp585_surface, temp_hindcast_surface) |>
  mutate(exposure_factor = "SST", factor_name = "Sea surface temperature")

BT <- create_anomaly_roms(temp_ssp585_bottom, temp_hindcast_bottom) |>
  mutate(exposure_factor = "BT", factor_name = "Bottom temperature")

SS <- create_anomaly_roms(salt_ssp585_surface, salt_hindcast_surface) |>
  mutate(exposure_factor = "SS", factor_name = "Surface salinity")

BS <- create_anomaly_roms(salt_ssp585_bottom, salt_hindcast_bottom) |>
  mutate(exposure_factor = "BS", factor_name = "Bottom salinity")

PhT <- create_anomaly_roms(PhT_ssp585_surface, PhT_hindcast_surface) |>
    mutate(exposure_factor = "PhT", factor_name = "Phytoplankton concentration")

ZP <- create_anomaly_roms(ZP_ssp585_surface, ZP_hindcast_surface) |>
    mutate(exposure_factor = "ZP", factor_name = "Zooplankton concentration")

# calculate anomaly for ESM variables (note different function)

SpH <- create_anomaly_gfdl("SpH", kriged[["pH_ssp585_surface"]], kriged[["pH_historical_surface"]]) |>
    mutate(exposure_factor = "SpH", factor_name = "Surface pH")

BpH <- create_anomaly_gfdl("BpH", kriged[["pH_ssp585_bottom"]], kriged[["pH_historical_bottom"]]) |>
    mutate(exposure_factor = "BpH", factor_name = "Bottom pH")

SO2 <- create_anomaly_gfdl("SO2", kriged[["o2_ssp585_surface"]], kriged[["o2_historical_surface"]]) |>
    mutate(exposure_factor = "SO2", factor_name = "Surface oxygen concentration")

BO2 <- create_anomaly_gfdl("BO2", kriged[["o2_ssp585_bottom"]], kriged[["o2_historical_bottom"]]) |>
    mutate(exposure_factor = "BO2", factor_name = "Bottom oxygen concentration")

AT <- create_anomaly_gfdl("AT", kriged[["AT_ssp585"]], kriged[["AT_historical"]]) |>
    mutate(exposure_factor = "AT", factor_name = "Air temperature")

PR <- create_anomaly_gfdl("PR", kriged[["PR_ssp585"]], kriged[["PR_historical"]]) |>
    mutate(exposure_factor = "PR", factor_name = "Precipitation")
  
# calculate anomaly for chosen MHW index (note different function)
MHW <- create_anomaly_mhw(mhw_ssp585, mhw_hindcast) |>
  mutate(exposure_factor = "MHW", factor_name = "Marine heatwaves")

# list of exposure factor data frames
exposure_factors_df <- list(SST, BT, SS, BS, PhT, ZP, SpH, BpH, SO2, BO2, AT, PR, MHW)
names(exposure_factors_df) <- c("SST", "BT", "SS", "BS", "PhT", "ZP", "SpH", "BpH", "SO2", "BO2", "AT", "PR", "MHW")
all_exposure <- rbind(SST, BT, SS, BS, PhT, ZP, SpH, BpH, SO2, BO2, AT, PR, MHW)

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

# create single anomaly and variability plots for all exposure factors and save locally 
for (k in 1:length(exposure_factors_df)){
  create_exposure_plot(exposure_factors_df[[k]], factor_full_names[k])
}

# create series of plots for anomalies: average mean, standard deviation, future change, and anomaly (13 plots)
for (k in 1:length(exposure_factors_df)){
  create_all_anomaly_plots(exposure_factors_df[[k]], exposure_factors[k])
}

# create series of plots for variability: future variance, historical variance, and variability (13 plots)
for (k in 1:length(exposure_factors_df)){
  create_all_variability_plots(exposure_factors_df[[k]], exposure_factors[k])
}

# ============================================================================
# SECTION 4: Species distributions
# ============================================================================
# Source and plot distribution maps for all species.
# ============================================================================

# file paths for species distribution data 
gdb_path <- here::here("data", "GOA_groundfish_2023.gdb") # NOAA Essential Fish Habitat maps 
scallop_path <- here::here("data", "EFH_2018_Scallop.gdb") # NOAA Essential Fish Habitat maps (presence-absence only)
bts_sdm_path <- here::here("data", "bts_sdms/") # estimates from the NOAA bottom trawl survey
diet_derived_path <- here::here("data", "diet_sdms/") # estimates from predator diet data
depth_temp_path <- here::here("data", "depth_temp_sdms/") # estimates based on species depth and temperature preferences

# list of all species, associated distribution file/layer, and exposure assignments
layer_names <- read_excel(here("exposure_assignments.xlsx"), sheet = "CVA_exposure_assignments", skip = 1,
                           col_names = c("group", "species_name", "scientific_name", "path", 
                                         "EFH_level", "layer", "sea_temperature", 
                                         "salinity", "ph", "phytoplankton", "zooplankton",
                                         "oxygen", "air_temperature", "precipitation", 
                                         "marine_heatwave", "data_quality")) |>
  drop_na(layer)

# clean layer_names and make lists
layer_names$air_temperature[layer_names$air_temperature == "NA"] <- NA
layer_names$precipitation[layer_names$precipitation == "NA"] <- NA
species_layers <- layer_names$layer
species_name <- layer_names$species_name
paths <- layer_names$path
EFH_level <- layer_names$EFH_level

# # test distribution plots
# plot_species_distribution("Sablefish")

# for loop to create plots for all species (57 plots)
for (i in 1:length(species_layers)) {
  this_species_name <- species_name[i]
  plot_species_distribution(this_species_name)
}

# ============================================================================
# SECTION 4: Calculate exposure scores
# ============================================================================
# Create overlap of species distribution and ROMS outputs. Calculate exposure 
# scores for anomaly and variability exposure factors. Save in results/ folder 
# for vulnerability calculation.
# ============================================================================

# create data frame to place exposure scores in 
exposure_scores <- layer_names |>
  dplyr::select(!c(path, EFH_level, layer)) |>
  pivot_longer(cols = c("sea_temperature", "salinity", "ph", "phytoplankton", "zooplankton",
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
  exposure_factors_list <- layer_names[i,]
  exposure_factors_list <- exposure_factors_list |> 
    dplyr::select(sea_temperature:marine_heatwave) # identify species in layer list 
  
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

# organize exposure scores
all_exposure_scores <- rbind(anomaly_scores, variability_scores) |>
  mutate(exposure_factor_full = recode(exposure_factor,
                                       "SST" = "Sea surface temperature",
                                       "BT" = "Bottom temperature",
                                       "SS" = "Surface salinity",
                                       "BS" = "Bottom salinity",
                                       "PhT" = "Phytoplankton concentration",
                                       "ZP" = "Zooplankton concentration",
                                       "SpH" = "Surface pH",
                                       "BpH" = "Bottom pH",
                                       "SO2" = "Surface oxygen concentration",
                                       "BO2" = "Bottom oxygen concentration",
                                       "AT" = "Air temperature",
                                       "PR" = "Precipitation",
                                       "MHW" = "Marine heatwaves"),
         name_full = paste(exposure_factor_full, type)) |>
  mutate(group = recode(group,
                        "Mollusc" = "Molluscs",
                        "Crustacean" = "Crustaceans",
                        "Elasmobranch" = "Elasmobranchs",
                        "Salmon" = "Salmonids",
                        "Forage" = "Forage fishes",
                        "Groundfish" = "Other groundfishes",
                        "Rockfish" = "Rockfishes",
                        "Flatfish" = "Flatfishes"))

all_exposure_scores$group  = ordered(all_exposure_scores$group,
                                 levels = c("Molluscs",
                                            "Crustaceans",
                                            "Elasmobranchs",
                                            "Salmonids",
                                            "Forage fishes",
                                            "Other groundfishes",
                                            "Rockfishes",
                                            "Flatfishes"))

all_exposure_scores <- all_exposure_scores |>
  arrange(group)

# save exposure results to csv file
write.csv(all_exposure_scores, "results/exposure_results.csv", row.names = FALSE)

# prep exposure data for bootstrapping with spread of tallies for exposure factor
n <- 1
bootstrap_anomalies <- exposure_scores
bootstrap_variability <- exposure_scores
for (i in 1:length(species_layers)) {
  # identify appropriate list of exposure factors for each species
  exposure_factors_list <- layer_names[i,]
  exposure_factors_list <- exposure_factors_list |> 
    dplyr::select(sea_temperature:marine_heatwave) # identify species in layer list 
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])
  
  for(k in 1:length(these_exposure_factors$V1)){
    original_exposure_data <- create_overlap(species_name[i], exposure_factors_df[[paste(these_exposure_factors$V1[k])]], these_exposure_factors$V1[k])
    
    anomaly_tallies <- assign_exposure_levels(original_exposure_data, "anomaly") |>
      dplyr::select(!c(total, prop)) |>
      pivot_wider(names_from = exposure_score,
                  values_from = count)
    variability_tallies <- assign_exposure_levels(original_exposure_data, "variability") |>
      dplyr::select(!c(total, prop)) |>
      pivot_wider(names_from = exposure_score,
                  values_from = count)

    bootstrap_anomalies[n,"sum_low"] <- anomaly_tallies$low
    bootstrap_anomalies[n,"sum_moderate"] <- anomaly_tallies$moderate
    bootstrap_anomalies[n,"sum_high"] <- anomaly_tallies$high
    bootstrap_anomalies[n,"sum_vh"] <- anomaly_tallies$very_high
    bootstrap_anomalies[n,"type"] <- "anomaly"
    
    bootstrap_variability[n,"sum_low"] <- variability_tallies$low
    bootstrap_variability[n,"sum_moderate"] <- variability_tallies$moderate
    bootstrap_variability[n,"sum_high"] <- variability_tallies$high
    bootstrap_variability[n,"sum_vh"] <- variability_tallies$very_high
    bootstrap_variability[n,"type"] <- "variability"
    
    n <- n + 1
  }
}

bootstrap_df <- rbind(bootstrap_anomalies, bootstrap_variability)
write.csv(bootstrap_df, "results/exposure_results_with_tallies.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: Summary plots
# ============================================================================
# Create summary plots of exposure scores for select or all species.
# ============================================================================

# plot distribution of scores for all exposure factors for each species
  # test plot
  # exposure_histogram_series("Walleye pollock", "anomaly")

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
  #     exposure_factors_list <- layer_names[i,]
  #   exposure_factors_list <- exposure_factors_list |> 
  #     dplyr::select(sea_temperature:marine_heatwave) # identify species in layer list 
  #   these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  #   these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])
  #   
  #   for(k in 1:length(these_exposure_factors$V1)){
  #     create_exposure_plots(species_name[i], exposure_factors_df[[paste(these_exposure_factors$V1[k])]], these_exposure_factors$V1[k])
  #   }
  # }
