# ============================================================================
# Exposure script
# ============================================================================
# Purpose: This script calculate exposure for all species for which EFH maps are available and
# a single exposure factor (surface temperature).

# Clear work space and free up memory
rm(list = ls())
gc()

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, here, terra, readxl)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
here::i_am("exposure_EFH.R")

# Load custom analysis functions by Albi
source(here("functions.R"))
source(here("exposure_functions.R"))

# ============================================================================
# SECTION 1: Load and Filter Data - Historical and Future
# ============================================================================
# The Arrow package allows "lazy evaluation" - you can filter the data before
# loading it into RAM, which is much faster than loading everything first.
# ============================================================================

load_data <- function(depth, variable_name) {
variable <- open_dataset(here(paste("data/processed/all_scenarios_bias_corrected_", variable_name, "_", depth, ".parquet", sep =""))) # load projections
future <- variable |> filter(run == "ssp585") |> # future projections
  filter(date > as.Date("2030-01-01")) |>  # Restrict dates to later than 2030
  filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
  collect()  # NOW load the filtered data into RAM
hindcast <- variable |> filter(run == "hindcast") |> # hindcast 
  filter(date > as.Date("1993-01-01")) |>  # Restrict dates to later than 1991
  filter(date < as.Date("2019-12-31")) |>  # Restrict to 2020
  collect()  # NOW load the filtered data into RAM

assign(paste(variable_name, "_", depth, "_future", sep = ""), future, envir = .GlobalEnv)
assign(paste(variable_name, "_", depth, "_hindcast", sep = ""), hindcast, envir = .GlobalEnv)
}

load_data("surface", "temp")
load_data("bottom", "temp")
load_data("surface", "salt")
load_data("bottom", "salt")

exposure_factors <- c("SST", "BT", "SS", "BS")

# ============================================================================
# SECTION 2: Calculate Anomaly
# ============================================================================
# Calculate standardized anomaly: (future mean - historical mean) / historical SD.
# ============================================================================

# function that calculates an anomaly for all exposure factors in this script
# currently set to summer months only; can adjust these
create_anomaly <- function(future_data, hindcast_data) {
  # calculate climate anomaly (future mean - historical mean / historical standard deviation)
  calculate_anomaly <- future_data |> 
    filter(month == "7" | month == "8" | month == "9") |> # can filter by winter or summer months by changing month numbers
    summarize(average_future = mean(value_dc), .by = c(cell_id, lon_rho, lat_rho)) |>
    left_join(
      hindcast_data |> 
        filter(month == "7" | month == "8" | month == "9") |> # filter in hindcast as well
        summarize(average_hist = mean(value), sd_hist = sd(value), .by = c(cell_id, lon_rho, lat_rho)), # for hindcast, value=value_dc
      by = join_by(lon_rho, lat_rho)
    ) |>
    mutate(anomaly = (average_future-average_hist)/sd_hist) |> # calculate anomaly 
    select(!c(average_future, average_hist, sd_hist))
  
  anomalies <- calculate_anomaly |> 
    st_as_sf(coords = c("lon_rho", "lat_rho"))
  st_crs(anomalies)= 4326
  anomalies <- anomalies |>
    mutate( # set scoring categories
      anomaly_bins = case_when(anomaly >= -5 & anomaly < -2 ~"very high",
                               anomaly >= -2 & anomaly < -1.5 ~"high",
                               anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                               anomaly >= -0.5 & anomaly < 0.5 ~"low",
                               anomaly >= 0.5 & anomaly < 1.5 ~"moderate",
                               anomaly >= 1.5 & anomaly < 2 ~"high",
                               anomaly >=2 & anomaly <= 5 ~"very high")
    )
  return(anomalies)
}

SST_anomaly <- create_anomaly(temp_surface_future, temp_surface_hindcast)
SST_anomaly <- SST_anomaly |> mutate(type = "SST")

BT_anomaly <- create_anomaly(temp_bottom_future, temp_bottom_hindcast)
BT_anomaly <- BT_anomaly |> mutate(type = "BT")

SS_anomaly <- create_anomaly(salt_surface_future, salt_surface_hindcast)
SS_anomaly <- SS_anomaly |> mutate(type = "SS")

BS_anomaly <- create_anomaly(salt_bottom_future, salt_bottom_hindcast)
BS_anomaly <- BS_anomaly |> mutate(type = "BS")

anomaly <- list(SST_anomaly, BT_anomaly, SS_anomaly, BS_anomaly) # combine all anomaly data frames

# Load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

# create anomaly plots for all exposure factors
# custom function from "exposure_functions.R"
for (k in 1:length(exposure_factors)){
  create_anomaly_plot(anomaly[[k]], exposure_factors[k])
}

# ============================================================================
# SECTION 3: Calculate exposure
# ============================================================================
# Map EFH maps and create overlap with ROMS outputs. Determine exposure and 
# create associated plots.
# ============================================================================

####### Load, Adjust, and Plot EFH Maps ######

#tell R where files are 
gdb_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/GOA_groundfish_2023.gdb"
# scallop_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/EFH_2018_Scallop.gdb"
# salmon_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/Salmon_2023.gdb"

#list all vector layers in the .gbd
layers <- vector_layers(gdb_path)
print(layers)

# scallop_layers <- vector_layers(scallop_path)
# print(scallop_layers)
# salmon_layers <- vector_layers(salmon_path)
# print(salmon_layers)

# get list of all the species layers we want to put through the function
layer_names <- read_excel(here("goa_efh_spp_lifestages.xlsx"), sheet = "groundfish",
                           col_names = c("species_name", "layer_name"))
species_layers<-unique(layer_names$layer_name)
species_name<-unique(layer_names$species_name)

# testing
plot_EFH_layer("GOA_adult_alaskaplaice_efh_level2_abundance_summer")

# for loop to create EFH plots for all species (30) 
# custom function from exposure_functions.R
for (i in 1:length(species_layers)) {
  plot_EFH_layer(species_layers[i])
}

####### Calculate exposure ######

# source custom functions from exposure_functions.R

# tests - ak plaice
create_exposure_plots("GOA_adult_alaskaplaice_efh_level2_abundance_summer", "Alaska plaice", SST_anomaly, "SST")
calculate_exposure_score("GOA_adult_alaskaplaice_efh_level2_abundance_summer", "Alaska plaice", SST_anomaly, "SST")

# create exposure plots for all EFH species and exposure factors #######
# create_exposure_plots is a custom function from exposure_functions.R
for (i in 1:length(species_layers)) {
  for(k in 1:length(exposure_factors)){
    create_exposure_plots(species_layers[i], species_name[i], anomaly[[k]], exposure_factors[k])
  }
}

# calculate exposure scores for all species and exposure factors and place them in the layer_names df
# calculate_exposure_score is a custom function from exposure_functions.R
for (i in 1:length(species_layers)) {
  for(k in 1:length(exposure_factors)){
  exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], anomaly[[k]], exposure_factors[k])
  layer_names[i,exposure_factors[k]] <- exposure_score
}
}

