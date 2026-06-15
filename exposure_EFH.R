# ============================================================================
# Exposure script
# ============================================================================
# Purpose: This script calculate exposure for all species for which EFH maps are available 
# and all exposure factors pulled from ROMS and GFDL ESM.

# Clear work space and free up memory
rm(list = ls())
gc()

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, here, terra, readxl, data.table)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
here::i_am("exposure_EFH.R")

# Load custom analysis functions by Albi
source(here("functions.R"))
source(here("exposure_functions.R"))

# ============================================================================
# SECTION 1A: Load and Filter ROMS Data - Historical and Future
# ============================================================================
# The Arrow package allows "lazy evaluation" - you can filter the data before
# loading it into RAM, which is much faster than loading everything first.
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

#combining phytoplankton (small and large)
PhT_ssp585_surface <- rbind(PhL_ssp585_surface, PhS_ssp585_surface) 
PhT_ssp585_surface <- PhT_ssp585_surface |> summarize(value_dc = sum(value_dc), .by = c(cell_id, lon_rho, lat_rho, year, month))
PhT_hindcast_surface <- rbind(PhL_hindcast_surface, PhS_hindcast_surface) 
PhT_hindcast_surface <- PhT_hindcast_surface |> summarize(value = sum(value), .by = c(cell_id, lon_rho, lat_rho, year, month))

#combining zooplankton (all)
ZP_ssp585_surface <- rbind(Cop_ssp585_surface, NCa_ssp585_surface, Eup_ssp585_surface, MZL_ssp585_surface, MZS_ssp585_surface) 
ZP_ssp585_surface <- PhT_ssp585_surface |> summarize(value_dc = sum(value_dc), .by = c(cell_id, lon_rho, lat_rho, year, month))
ZP_hindcast_surface <- rbind(Cop_hindcast_surface, NCa_hindcast_surface, Eup_hindcast_surface, MZL_hindcast_surface, MZS_hindcast_surface) 
ZP_hindcast_surface <- PhT_hindcast_surface |> summarize(value = sum(value), .by = c(cell_id, lon_rho, lat_rho, year, month))

# ============================================================================
# SECTION 1B: Load GFDL Data - Historical and Future
# ============================================================================
# Load parquet files from data folder.
# These files are created using the load_gfdl_data.R script.
# ============================================================================

ph_ssp585 <- open_dataset(here("data/pH/ph_ssp585_surface.parquet"))
ph_historical <- open_dataset(here("data/pH/ph_historical_surface_2020.parquet"))

o2_ssp585 <- open_dataset(here("data/o2/o2_ssp585_surface.parquet"))
o2_historical <- open_dataset(here("data/o2/o2_historical_surface_2020.parquet"))

at_ssp585 <- open_dataset(here("data/tas/airtemp_ssp585_na.parquet"))
at_historical <- open_dataset(here("data/tas/airtemp_historical_na_2020.parquet"))

pr_ssp585 <- open_dataset(here("data/pr/precip_ssp585_na.parquet"))
pr_historical <- open_dataset(here("data/pr/precip_historical_na_2020.parquet"))

exposure_factors <- c("SST", "BT", "SS", "BS", "PhT", "ZP", "PH", "O2", "AT", "PR") 
# sea-surface temp, bottom temp, surface salinity, bottom salinity,
# phytoplankton concentration, zooplankton concentration, pH (ocean acidification),
# oxygen concentration, air temperature, precipitation

# ============================================================================
# SECTION 2: Calculate Anomalies
# ============================================================================
# Calculate standardized anomalies: (future mean - historical mean) / historical SD
# for each exposure factor. Functions from exposure_functions.R.
# ============================================================================

# calculates an anomaly for all exposure factors in this script
# currently set to summer months only; can adjust these in exposure_functions.R

SST_anomaly <- create_anomaly_roms(temp_ssp585_surface, temp_hindcast_surface) # sea surface temp
SST_anomaly <- SST_anomaly |> mutate(type = "SST")

BT_anomaly <- create_anomaly_roms(temp_ssp585_bottom, temp_hindcast_bottom) # bottom temp
BT_anomaly <- BT_anomaly |> mutate(type = "BT")

SS_anomaly <- create_anomaly_roms(salt_ssp585_surface, salt_hindcast_surface) # surface salinity
SS_anomaly <- SS_anomaly |> mutate(type = "SS")

BS_anomaly <- create_anomaly_roms(salt_ssp585_bottom, salt_hindcast_bottom) # bottom salinity
BS_anomaly <- BS_anomaly |> mutate(type = "BS")

PhT_anomaly <- create_anomaly_roms(PhT_ssp585_surface, PhT_hindcast_surface) # phytoplankton
PhT_anomaly <- PhT_anomaly |> mutate(type = "PhT")

ZP_anomaly <- create_anomaly_roms(ZP_ssp585_surface, ZP_hindcast_surface) # zooplankton
ZP_anomaly <- ZP_anomaly |> mutate(type = "ZP")

# note different function for GFDL exposure factors
PH_anomaly <- create_anomaly_gfdl("PH", ph_ssp585, ph_historical) # pH
PH_anomaly <- PH_anomaly |> mutate(type = "PH")

O2_anomaly <- create_anomaly_gfdl("O2", o2_ssp585, o2_historical) # oxygen
O2_anomaly <- O2_anomaly |> mutate(type = "O2")

AT_anomaly <- create_anomaly_gfdl("AT", at_ssp585, at_historical) # air temperature
AT_anomaly <- AT_anomaly |> mutate(type = "AT")

PR_anomaly <- create_anomaly_gfdl("PR", pr_ssp585, pr_historical) # precipitation
PR_anomaly <- PR_anomaly |> mutate(type = "PR")

anomaly <- list(SST_anomaly, BT_anomaly, SS_anomaly, BS_anomaly, PhT_anomaly, ZP_anomaly, PH_anomaly, O2_anomaly, AT_anomaly, PR_anomaly) # combine all anomaly data frames
names(anomaly) <- c("SST_anomaly", "BT_anomaly", "SS_anomaly", "BS_anomaly", "PhT_anomaly", "ZP_anomaly", "PH_anomaly", "O2_anomaly", "AT_anomaly", "PR_anomaly")

sf::sf_use_s2(FALSE) # IMPORTANT - turns off spherical geometry

# Load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

# load GOA shading
GOA = ne_countries(scale = "medium", 
                     returnclass = "sf") |>
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude 

# load Canada shading
canada = ne_countries(scale = "medium", country = "Canada", returnclass = "sf") |>
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude 

create_anomaly_plot(AT_anomaly, "AT")

# create anomaly plots for all exposure factors and save locally
# function from exposure_functions.R
for (k in 1:length(exposure_factors)){
  create_anomaly_plot(anomaly[[k]], exposure_factors[k])
}

# ============================================================================
# SECTION 3: Load and plot EFH maps.
# ============================================================================
# Source and create EFH maps for each species by transforming .gdb files into 
# sf objects. 
# ============================================================================

# tell R where EFH files are 
gdb_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/GOA_groundfish_2023.gdb"
scallop_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/EFH_2018_Scallop.gdb"
bts_sdm_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/sdms/"
# salmon_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/Salmon_2023.gdb"

# list all vector layers in the .gdb
layers <- vector_layers(gdb_path)
print(layers)

# list of all the species layers we have EFH for 
layer_names <- read_excel(here("goa_efh_spp_lifestages.xlsx"), sheet = "groundfish", skip = 1,
                           col_names = c("group", "species_name", "path", "EFH_level", "layer", "sea_temperature", 
                                         "salinity", "ph", "phytoplankton", "zooplankton",
                                         "oxygen", "air_temperature", "precipitation")) |>
  drop_na(layer)
layer_names$air_temperature[layer_names$air_temperature == "NA"] <- NA
layer_names$precipitation[layer_names$precipitation == "NA"] <- NA
species_layers <- unique(layer_names$layer)
species_name <- unique(layer_names$species_name)
paths <- layer_names$path
EFH_level <- layer_names$EFH_level

# test plot
plot_EFH_layer("gdb_path", "GOA_adult_alaskaplaice_efh_level2_abundance_summer", "2")
plot_EFH_layer("scallop_path", "_Weathervane_scallop_adult_EFH_Level1", "1")
plot_EFH_layer("bts_sdm_path", "predictions_Pacific.capelin.rda", "2")

# for loop to create EFH plots for all species (30) 
# function from exposure_functions.R
for (i in 1:length(species_layers)) {
  this_path_name <- paths[i]
  this_species_layer <- species_layers[i]
  this_EFH_level <- EFH_level[i]
  plot_EFH_layer(this_path_name, this_species_layer, this_EFH_level)
}

# ============================================================================
# SECTION 4: Calculate exposure.
# ============================================================================
# Create overlap of EFH and ROMS outputs. Determine exposure and 
# create associated plots.
# ============================================================================

# functions from "exposure_functions.R" script
# create_exposure_plots() creates 3 plots and saves locally
# 1. exposure map (binned anomalies within species distribution)
# 2. histogram of anomalies across the species distribution
# 3. distribution of exposure scores and weighted average 
# calculate_exposure_score() calculates exposure scores and puts them in the layer_names data frame

capelin <- create_overlap("bts_sdm_path", "predictions_Pacific.capelin.rda", "2", "Capelin", SST_anomaly, "SST")
alaska_plaice <- create_overlap("gdb_path", "GOA_adult_alaskaplaice_efh_level2_abundance_summer", "2", "Alaska plaice", BT_anomaly, "BT")

# test out some plots
create_exposure_plots("gdb_path", "GOA_adult_alaskaplaice_efh_level2_abundance_summer", "2", "Alaska plaice", BT_anomaly, "BT")
create_exposure_plots("GOA_adult_alaskaskate_efh_level2_abundance_summer", "Alaska skate", PH_anomaly, "PH")
create_exposure_plots("GOA_adult_starryflounder_efh_level2_abundance_summer", "Starry flounder", AT_anomaly, "AT")
create_exposure_plots("bts_sdm_path", "predictions_Pacific.capelin.rda", "2", "Capelin", PH_anomaly, "PH")

# example code
 exposure_factors_list <- layer_names[20,6:13]
 these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
 these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])

# create 3 exposure plots for all EFH species and assigned exposure factors
for (i in 1:length(species_layers)) {
  # identify appropriate list of exposure factors for each species
  exposure_factors_list <- layer_names[i,6:13]
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])
  
  this_path_name <- paths[i]
  this_species_layer <- species_layers[i]
  this_EFH_level <- EFH_level[i]
  for(k in 1:length(these_exposure_factors$V1)){
    create_exposure_plots(this_path_name, this_species_layer, this_EFH_level, species_name[i], anomaly[[paste(these_exposure_factors$V1[k], "_anomaly", sep = "")]], these_exposure_factors$V1[k])
  }
}

# exposure_scores <- data.frame(species = character(), exposure_factor = character(), score = numeric())

# create data frame
exposure_scores <- layer_names |>
  dplyr::select(!c(path, EFH_level, layer)) |>
  pivot_longer(cols=c("sea_temperature", "salinity", "ph", "phytoplankton", "zooplankton",
                      "oxygen", "air_temperature", "precipitation"),
                               names_to = "variable",
                               values_to = "exposure_factor") |>
  dplyr::select(!variable) |>
  drop_na(exposure_factor) 

# calculate exposure scores for all species and exposure factors and place them in the df
row_pointer <- 1

for (i in 1:length(species_layers)) {
  # identify appropriate list of exposure factors for each species
  exposure_factors_list <- layer_names[i,6:13]
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])

  for(k in 1:length(these_exposure_factors$V1)){
    original_exposure_data <- create_overlap(paths[i], species_layers[i], EFH_level[i], species_name[i], anomaly[[paste(these_exposure_factors$V1[k], "_anomaly", sep = "")]], these_exposure_factors$V1[k])
    score <- calculate_exposure_score(original_exposure_data, paths[i], species_layers[i], EFH_level[i], species_name[i], anomaly[[paste(these_exposure_factors$V1[k], "_anomaly", sep = "")]], these_exposure_factors$V1[k])
  
    exposure_scores[row_pointer,"score"] <- score
  
    row_pointer <- row_pointer + 1
  }
}

# save csv with exposure scores
write.csv(exposure_scores, "exposure_factor_scores.csv", row.names = FALSE)

# format exposure scores for vulnerability calculation below 

# # reformat data frame
# indv_exposure_results <- layer_names |>
#   pivot_longer(cols=c("SST", "BT", "SS", "BS", "PhT", "ZP", "PH", "O2", "AT", "PR"),
#                names_to = "exposure_factor_short",
#                values_to = "score") |>
#   mutate(exposure_factor = recode(exposure_factor_short,
#                                   "SST" = "Sea surface temperature",
#                                   "BT" = "Bottom temperature",
#                                   "SS" = "Surface salinity",
#                                   "BS" = "Bottom salinity",
#                                   "PhT" = "Phytoplankton concentration",
#                                   "ZP" = "Zooplankton concentration",
#                                   "PH" = "pH",
#                                   "O2" = "Oxygen concentration",
#                                   "AT" = "Air temperature",
#                                   "PR" = "Precipitation"))
# 
# # apply logic rule
# exposure_results <- indv_exposure_results |>
#   group_by(group, species_name) |>
#   summarize(
#     above_3.5 = sum(score >= 3.5), # number of mean scores above 3.5 - "very high" 
#     above_3 = sum(score >= 3), # number of mean scores above 3 - "high" 
#     above_2.5 = sum(score >=2.5), # number of mean scores above 2.5 - "moderate" 
#     exposure = ifelse(above_3.5 >= 3, "very high", ifelse(above_3 >= 2, "high", ifelse(above_2.5 >= 2, "moderate", "low"))),
#   ) |>
#   ungroup() |>
#   dplyr::select(!c(above_3.5, above_3, above_2.5))
# 
# # apply ALTERNATIVE logic rule
# alt_exposure_results <- indv_exposure_results |>
#   group_by(group, species_name) |>
#   summarize(
#     above_3.5 = sum(score >= 3.5),
#     above_2.5 = sum(score >= 2.5), # mean scores between 2.5 and 3.5 would be considered as “high"
#     above_1.5 = sum(score >=1.5), # mean scores between 1.5 and 2.5 would be considered as “moderate”
#     exposure = ifelse(above_3.5 >= 3, "very high", ifelse(above_2.5 >= 2, "high", ifelse(above_1.5 >= 2, "moderate", "low"))),
#   ) |>
#   ungroup() |>
#   dplyr::select(!c(above_3.5, above_2.5, above_1.5))
# 
# exposure_results$exposure  = ordered(exposure_results$exposure, 
#                                      levels = c("low",
#                                                 "moderate",
#                                                 "high",
#                                                 "very high"))
# # add number for each exposure score
# exposure_results <- exposure_results |>
#   mutate(exposure_number = recode(exposure, 
#                                   "low" = "1",
#                                   "moderate" = "2", 
#                                   "high" = "3", 
#                                   "very high" = "4"))
# 
# # save csv of final exposure scores for vulnerability calculation
# write.csv(exposure_results, "exposure_scores_final.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: More plots
# ============================================================================
# Create plot that includes distribution of exposure scores for all exposure 
# factors for each species.
# ============================================================================

# plot distribution of exposure scores for all exposure factors for each species
sp_exposure_histograms <- function(path, species_layer, EFH_level, species_name){
  exposure_plot_list = list() # create list
  # select correct exposure factors for each species
  exposure_factors_list <- layer_names[i,6:13]
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])
  
  for(i in 1:length(these_exposure_factors$V1)){
    original_exposure_data <- create_overlap(path, species_layer, EFH_level, species_name, anomaly[[paste(these_exposure_factors$V1[i], "_anomaly", sep = "")]], these_exposure_factors$V1[i])
    exposure_plot <- assign_exposure_levels(original_exposure_data)
    exposure_plot <- exposure_plot |>
      mutate(exposure_factor = these_exposure_factors$V1[i])
    # assign(paste(exposure_factors[i], "_exposure_plot", sep = ""), exposure_plot, envir = .GlobalEnv)
    exposure_plot_list[[i]] <- exposure_plot # add it to the list
  }
  # create single data frame with every df from the for loop
  # combine all exposure factor results into a single df per species
  group_exposure_plot <- do.call(rbind, exposure_plot_list) 
  
  # exp_fact_mean <- calculate_exposure_score(species_layer, species_name, anomaly_data, exposure_factor_name)
  # find way to add specific exposure scores to each facet
  
  group_plot <- ggplot(group_exposure_plot) +
    geom_col(mapping = aes(x = exposure_score, y = prop, fill = exposure_score), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
    facet_wrap(~exposure_factor) +
    scale_x_discrete(labels = c("L", "M", "H", "V")) +
    scale_y_continuous(labels = scales::percent) +
    ylab("Percent") +
    xlab("Exposure Score") +
    scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
    #annotate("text", x = I(0.8), y = I(0.8), label = paste("Exposure = ", exp_fact_mean, sep = "")) + # add exposure score in top right corner
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  
  ggsave(paste(species_name, "all_exposure_scores.png", sep="_"), path = here("plots/Exposure plots/"), plot = group_plot, device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)
}

# test plot
sp_exposure_histograms("gdb_path", "GOA_adult_walleyepollock_efh_level2_abundance_summer", "2", "Walleye pollock")

# for to create plot for all species
for (i in 1:length(species_layers)) {
  sp_exposure_histograms(paths[i], species_layers[i], EFH_level[i], species_name[i])
}
