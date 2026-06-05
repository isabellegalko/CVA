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

# tell R where EFH  files are 
gdb_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/GOA_groundfish_2023.gdb"
# scallop_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/EFH_2018_Scallop.gdb"
# salmon_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/Salmon_2023.gdb"

# list all vector layers in the .gdb
layers <- vector_layers(gdb_path)
print(layers)

# list of all the species layers we have EFH for 
layer_names <- read_excel(here("goa_efh_spp_lifestages.xlsx"), sheet = "groundfish",
                           col_names = c("group", "species_name", "layer_name"))
species_layers <- unique(layer_names$layer_name)
species_name <- unique(layer_names$species_name)

# test plot
plot_EFH_layer("GOA_adult_alaskaplaice_efh_level2_abundance_summer")

# for loop to create EFH plots for all species (30) 
# function from exposure_functions.R
for (i in 1:length(species_layers)) {
  plot_EFH_layer(species_layers[i])
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

# test out some plots
create_exposure_plots("GOA_adult_alaskaplaice_efh_level2_abundance_summer", "Alaska plaice", SST_anomaly, "SST")
create_exposure_plots("GOA_adult_alaskaskate_efh_level2_abundance_summer", "Alaska skate", PH_anomaly, "PH")
create_exposure_plots("GOA_adult_starryflounder_efh_level2_abundance_summer", "Starry flounder", AT_anomaly, "AT")

# create 3 exposure plots for all EFH species and exposure factors
for (i in 1:length(species_layers)) {
  for(k in 1:length(exposure_factors)){
    create_exposure_plots(species_layers[i], species_name[i], anomaly[[k]], exposure_factors[k])
  }
}

# calculate exposure scores for all species and exposure factors and place them in the layer_names df
for (i in 1:length(species_layers)) {
  for(k in 1:length(exposure_factors)){
  exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], anomaly[[k]], exposure_factors[k])
  layer_names[i,exposure_factors[k]] <- exposure_score
}
}

# save csv with exposure scores
write.csv(layer_names, "exposure_factor_scores.csv", row.names = FALSE)

# reformat data frame
indv_exposure_results <- layer_names |>
  pivot_longer(cols=c("SST", "BT", "SS", "BS", "PhT", "ZP", "PH", "O2", "AT", "PR"),
               names_to = "exposure_factor_short",
               values_to = "score") |>
  mutate(exposure_factor = recode(exposure_factor_short,
                                  "SST" = "Sea surface temperature",
                                  "BT" = "Bottom temperature",
                                  "SS" = "Surface salinity",
                                  "BS" = "Bottom salinity",
                                  "PhT" = "Phytoplankton concentration",
                                  "ZP" = "Zooplankton concentration",
                                  "PH" = "pH",
                                  "O2" = "Oxygen concentration",
                                  "AT" = "Air temperature",
                                  "PR" = "Precipitation"))

# apply logic rule
exposure_results <- indv_exposure_results |>
  group_by(group, species_name) |>
  summarize(
    above_3.5 = sum(score >= 3.5), # number of mean scores above 3.5 - "very high" 
    above_3 = sum(score >= 3), # number of mean scores above 3 - "high" 
    above_2.5 = sum(score >=2.5), # number of mean scores above 2.5 - "moderate" 
    exposure = ifelse(above_3.5 >= 3, "very high", ifelse(above_3 >= 2, "high", ifelse(above_2.5 >= 2, "moderate", "low"))),
  ) |>
  ungroup() |>
  select(!c(above_3.5, above_3, above_2.5))

# apply ALTERNATIVE logic rule
alt_exposure_results <- indv_exposure_results |>
  group_by(group, species_name) |>
  summarize(
    above_3.5 = sum(score >= 3.5),
    above_2.5 = sum(score >= 2.5), # mean scores between 2.5 and 3.5 would be considered as “high"
    above_1.5 = sum(score >=1.5), # mean scores between 1.5 and 2.5 would be considered as “moderate”
    exposure = ifelse(above_3.5 >= 3, "very high", ifelse(above_2.5 >= 2, "high", ifelse(above_1.5 >= 2, "moderate", "low"))),
  ) |>
  ungroup() |>
  select(!c(above_3.5, above_2.5, above_1.5))

exposure_results$exposure  = ordered(exposure_results$exposure, 
                                     levels = c("low",
                                                "moderate",
                                                "high",
                                                "very high"))
# add number for each exposure score
exposure_results <- exposure_results |>
  mutate(exposure_number = recode(exposure, 
                                  "low" = "1",
                                  "moderate" = "2", 
                                  "high" = "3", 
                                  "very high" = "4"))

# save csv of final exposure scores for vulnerability calculation
write.csv(exposure_results, "exposure_scores_final.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: More plots
# ============================================================================
# 
# ============================================================================

# exposure_plot <- assign_exposure_levels("GOA_adult_walleyepollock_efh_level2_abundance_summer", "Walleye pollock", SST_anomaly, "SST")

# plot distribution of exposure scores for all exposure factors for each species
sp_exposure_histograms <- function(species_layer, species_name){
  exposure_plot_list = list() # create list
  for(i in 1:length(exposure_factors)){
    exposure_plot <- assign_exposure_levels(species_layer, species_name, anomaly[[i]], exposure_factors[i])
    exposure_plot <- exposure_plot |>
      mutate(exposure_factor = exposure_factors[i])
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
  
  ggsave(paste(species_layer, "exposure_scores.png", sep="_"), path = here(paste("plots/Exposure plots/", species_name, sep ="")), plot = group_plot, device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)
}

# test
sp_exposure_histograms("GOA_adult_walleyepollock_efh_level2_abundance_summer", "Walleye pollock")

for (i in 1:length(species_layers)) {
  sp_exposure_histograms(species_layers[i], species_name[i])
}
