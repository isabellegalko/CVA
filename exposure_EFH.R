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
  pacman::p_load(ggplot2, tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, terra, readxl, data.table, gstat, ggpubr, cowplot, patchwork)
    pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
here::i_am("exposure_EFH.R")
source(here("exposure_functions.R")) # Load custom analysis functions by Isabelle

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
ZP_ssp585_surface <- ZP_ssp585_surface |> summarize(value_dc = sum(value_dc), .by = c(cell_id, lon_rho, lat_rho, year, month))
ZP_hindcast_surface <- rbind(Cop_hindcast_surface, NCa_hindcast_surface, Eup_hindcast_surface, MZL_hindcast_surface, MZS_hindcast_surface) 
ZP_hindcast_surface <- ZP_hindcast_surface |> summarize(value = sum(value), .by = c(cell_id, lon_rho, lat_rho, year, month))

# ============================================================================
# SECTION 1B: Load GFDL Data - Historical and Future
# ============================================================================
# Load parquet files from data folder. Plot variogram. Krige.
# These files are created using the load_gfdl_data.R script.
# ============================================================================

load("10km_grid.rda")
prediction.grid = as.data.frame(grid)
grid_sf = st_as_sf(prediction.grid, coords = c("lon_rho", "lat_rho"), crs = 4326) |>
  st_shift_longitude()

# fit and print variogram model results
# identify psill, range, and nugget
find_vgm_values <- function(vgm_data){
  # vgm_data = vgm_data |>
  #   summarize(mean_summer = mean(value), 
  #             .by = c(cell_id, lon, lat)) |> 
  #   sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  this_variable.vgm = variogram(mean_summer ~ 1, data = vgm_data)
  plot(this_variable.vgm)
  this_variable.vgm_fit = fit.variogram(this_variable.vgm, model=vgm(psill=0.1, model="Gau", range=450, nugget=0.1))
  plot(this_variable.vgm, this_variable.vgm_fit)
  
  # print and return table with values 
  print(this_variable.vgm_fit)
  return(this_variable.vgm_fit)
}

# kriges by year
# identifies the best fit psill, nugget, and range values for each year 
run_krige <- function(data){
  for_kriging <- data |>
    filter(month == c("7", "8", "9")) |> # filter by summer months 
    summarize(mean_summer = mean(value), .by = c(cell_id, lon, lat, year)) |> # need only one point per location to krige by year 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  kriged_list = list()
  for(i in unique(for_kriging$year)) {
    df = subset(for_kriging, year == i)
    
    # identify year-specific psill, nugget, and range values
    this_variable.vgm_fit <- find_vgm_values(df)
    
    nugget_val <- this_variable.vgm_fit[1, "psill"] # nugget
    sill_val   <- this_variable.vgm_fit[2, "psill"] # partial sill
    range_val  <- this_variable.vgm_fit[2, "range"] # range
    
    vgm = variogram(mean_summer ~ 1, data = df)
    vgm_fit = fit.variogram(vgm, model = vgm(psill = sill_val, model = "Gau", range = range_val, nugget = nugget_val))
    
    krig = gstat::krige(mean_summer ~ 1, df, grid_sf, model = vgm_fit)
    
    krig = as.data.frame(krig)
    krig$year = i
    kriged_list[[i]] = krig 
  }
  all_kriged = dplyr::bind_rows(kriged_list) |>
    rename(kriged = var1.pred)
  
  return(all_kriged)
}

pH_ssp585_surface <- open_dataset(here("data/pH/ph_ssp585_surface.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

pH_historical_surface <- open_dataset(here("data/pH/ph_historical_surface.parquet")) |>
  filter(year >= 2005 & year <= 2020) |>
    collect()

pH_ssp585_bottom <- open_dataset(here("data/pH/ph_ssp585_bottom.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

pH_historical_bottom <- open_dataset(here("data/pH/ph_historical_bottom.parquet")) |>
  filter(year >= 2005 & year <= 2020) |>
    collect()

o2_ssp585_surface <- open_dataset(here("data/o2/o2_ssp585_surface.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

o2_historical_surface <- open_dataset(here("data/o2/o2_historical_surface.parquet")) |>
  filter(year >= 2005 & year <= 2020) |>
    collect()

o2_ssp585_bottom <- open_dataset(here("data/o2/o2_ssp585_bottom.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

o2_historical_bottom <- open_dataset(here("data/o2/o2_historical_bottom.parquet")) |>
  filter(year >= 2005 & year <= 2020) |>
    collect()

AT_ssp585 <- open_dataset(here("data/tas/airtemp_ssp585_na.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

AT_historical <- open_dataset(here("data/tas/airtemp_historical_na_2020.parquet")) |>
  filter(year >= 2005 & year <= 2020) |>
    collect()

PR_ssp585 <- open_dataset(here("data/pr/precip_ssp585_na.parquet")) |>
  filter(year >= 2030 & year <= 2059) |>
    collect()

PR_historical <- open_dataset(here("data/pr/precip_historical_na_2020.parquet")) |>
  filter(year >= 2005 & year <= 2020) |>
    collect()

## NOTE: THERE ARE STILL MODEL CONVERGENCE ERRORS FOR PRECIPITATION!!!

gfdl_files <- c("pH_ssp585_surface", "pH_historical_surface", "pH_ssp585_bottom", "pH_historical_bottom", "o2_ssp585_surface", "o2_historical_surface", "o2_ssp585_bottom", "o2_historical_bottom", "AT_ssp585", "AT_historical", "PR_ssp585", "PR_historical")
  
kriged <- c()
for(i in 1:length(gfdl_files)){
  this_data = get(gfdl_files[i]) 
  kriged[[gfdl_files[i]]] <- run_krige(this_data)
}

exposure_factors <- c("SST", "BT", "SS", "BS", "PhT", "ZP", "SPH", "BPH", "SO2", "BO2", "AT", "PR") 
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

# calculate anomaly for ESM variables (note different function in exposure_functions.R)

SPH_anomaly <- create_anomaly_gfdl("SPH", kriged[["pH_ssp585_surface"]], kriged[["pH_historical_surface"]]) # surface pH
  SPH_anomaly <- SPH_anomaly |> mutate(type = "SPH")

BPH_anomaly <- create_anomaly_gfdl("BPH", kriged[["pH_ssp585_bottom"]], kriged[["pH_historical_bottom"]]) # bottom pH
  BPH_anomaly <- BPH_anomaly |> mutate(type = "BPH")

SO2_anomaly <- create_anomaly_gfdl("SO2", kriged[["o2_ssp585_surface"]], kriged[["o2_historical_surface"]]) # surface oxygen
  SO2_anomaly <- SO2_anomaly |> mutate(type = "SO2")

BO2_anomaly <- create_anomaly_gfdl("BO2", kriged[["o2_ssp585_bottom"]], kriged[["o2_historical_bottom"]]) # bottom oxygen
  BO2_anomaly <- BO2_anomaly |> mutate(type = "BO2")

AT_anomaly <- create_anomaly_gfdl("AT", kriged[["AT_ssp585"]], kriged[["AT_historical"]]) # air temperature
  AT_anomaly <- AT_anomaly |> mutate(type = "AT")

PR_anomaly <- create_anomaly_gfdl("PR", kriged[["PR_ssp585"]], kriged[["PR_historical"]]) # precipitation
  PR_anomaly <- PR_anomaly |> mutate(type = "PR")

anomaly <- list(SST_anomaly, BT_anomaly, SS_anomaly, BS_anomaly, PhT_anomaly, ZP_anomaly, SPH_anomaly, BPH_anomaly, SO2_anomaly, BO2_anomaly, AT_anomaly, PR_anomaly) # combine all anomaly data frames
names(anomaly) <- c("SST_anomaly", "BT_anomaly", "SS_anomaly", "BS_anomaly", "PhT_anomaly", "ZP_anomaly", "SPH_anomaly", "BPH_anomaly", "SO2_anomaly", "BO2_anomaly", "AT_anomaly", "PR_anomaly")

sf::sf_use_s2(FALSE) # IMPORTANT - turns off spherical geometry

# load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") |>
  st_crop(xmin = -172, xmax = -130, ymin = 50, ymax = 62) |>  # crop to GOA region
  st_shift_longitude()  # convert to 0-360° longitude to match ROMS data

# load GOA shading
GOA = ne_countries(scale = "medium", 
                   returnclass = "sf") |>
  st_crop(xmin = -172, xmax = -130, ymin = 50, ymax = 62) |>  # crop to GOA region
  st_shift_longitude()  # convert to 0-360° longitude 

# load Canada shading
canada = ne_countries(scale = "medium", country = "Canada", returnclass = "sf") |>
  st_crop(xmin = -172, xmax = -130, ymin = 50, ymax = 62) |>  # crop to GOA region
  st_shift_longitude()  # convert to 0-360° longitude 

# create single anomaly plots for all exposure factors and save locally
# function from exposure_functions.R
for (k in 1:length(exposure_factors)){
  create_anomaly_plot(anomaly[[k]], exposure_factors[k])
}

# create series of plots: show average mean, sd, future change, and anomalies
for (k in 1:length(exposure_factors)){
  create_all_anomaly_plots(anomaly[[k]], exposure_factors[k])
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
bts_sdm_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/bts_sdms/"
diet_derived_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/diet_sdms/"
depth_temp_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/depth_temp_sdms/"

# list all vector layers in the .gdb
layers <- vector_layers(gdb_path)
print(layers)

# list of all the species layers we have EFH for 
layer_names <- read_excel(here("goa_efh_spp_lifestages.xlsx"), sheet = "groundfish", skip = 1,
                           col_names = c("group", "species_name", "path", "EFH_level", "layer", "sea_temperature", 
                                         "salinity", "ph", "phytoplankton", "zooplankton",
                                         "oxygen", "air_temperature", "precipitation", "data_quality")) |>
  drop_na(layer)
layer_names$air_temperature[layer_names$air_temperature == "NA"] <- NA
layer_names$precipitation[layer_names$precipitation == "NA"] <- NA
species_layers <- layer_names$layer
species_name <- layer_names$species_name
paths <- layer_names$path
EFH_level <- layer_names$EFH_level

# # test plots
# plot_species_distribution("gdb_path", "GOA_adult_redbandedrockfish_efh_level2_abundance_summer", "2", "Redbanded rockfish")
# plot_species_distribution("scallop_path", "_Weathervane_scallop_adult_EFH_Level1", "1", "Weathervane scallop")
# plot_species_distribution("bts_sdm_path", "predictions_king.crabs.rda", "2", "Red king crab")
# plot_species_distribution("diet_derived_path", "predictions_Pandalid shrimps_diet_derived.rda", "2", "Spot shrimp")
# plot_species_distribution("depth_temp_path", "geoduck_clam_predictions.rda", "1", "Geoduck clam")

# for loop to create EFH plots for all species 
# function from exposure_functions.R
for (i in 1:length(species_layers)) {
  this_path_name <- paths[i]
  this_species_layer <- species_layers[i]
  this_EFH_level <- EFH_level[i]
  this_species_name <- species_name[i]
  plot_species_distribution(this_path_name, this_species_layer, this_EFH_level, this_species_name)
}

# ============================================================================
# SECTION 4: Calculate exposure.
# ============================================================================
# Create overlap of EFH and ROMS outputs. Determine exposure scores.
# ============================================================================

# test overlap function
# geoduck_clam <- create_overlap("depth_temp_path", "geoduck_clam_predictions.rda", "1", "Geoduck clam", BT_anomaly, "BT")

# create data frame to place exposure scores in 
exposure_scores <- layer_names |>
  dplyr::select(!c(path, EFH_level, layer)) |>
  pivot_longer(cols=c("sea_temperature", "salinity", "ph", "phytoplankton", "zooplankton",
                      "oxygen", "air_temperature", "precipitation"),
                               names_to = "variable",
                               values_to = "exposure_factor") |>
  dplyr::select(!variable) |>
    drop_na(exposure_factor) 

# test calculating exposure scores
# alaska_skate <- create_overlap("gdb_path", "GOA_adult_alaskaskate_efh_level2_abundance_summer", "2", "Alaska skate", BPH_anomaly, "BPH")
# calculate_exposure_score(alaska_skate, "gdb_path", "GOA_adult_alaskaskate_efh_level2_abundance_summer", "2", "Alaska skate", BPH_anomaly, "BPH")

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
# note: currently set to reference period of 2005-2020
write.csv(exposure_scores, "results/exposure_factor_scores_2005-2020.csv", row.names = FALSE)

# ============================================================================
# SECTION 5: Summary plots
# ============================================================================
# For each species, create plot that includes distribution of exposure scores 
# for all exposure factors for each species. 
# ============================================================================

# test plots
sp_exposure_histograms("gdb_path", "GOA_adult_walleyepollock_efh_level2_abundance_summer", "2", "Walleye pollock")
sp_exposure_histograms("depth_temp_path", "geoduck_clam_predictions.rda", "1", "Geoduck clam")

# create series of plots for all species
for (i in 1:length(species_layers)) {
  exposure_histogram_series(paths[i], species_layers[i], EFH_level[i], species_name[i])
}

create_exposure_plots("gdb_path", "GOA_adult_walleyepollock_efh_level2_abundance_summer", "2", "Walleye pollock", BT_anomaly, "BT")

# functions from "exposure_functions.R" script
# create_exposure_plots() creates 3 plots and saves locally
# 1. exposure map (binned anomalies within species distribution)
# 2. histogram of anomalies across the species distribution
# 3. distribution of exposure scores and weighted average 
# calculate_exposure_score() calculates exposure scores and puts them in the layer_names data frame
# create 3 exposure plots for all EFH species and assigned exposure factors (warning: creates ~900 plots)
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

# ============================================================================
# SECTION 6: Manuscript plots
# ============================================================================
# Create a plot that shows all exposure maps and histograms for a single species 
# as an example of the exposure process.
# ============================================================================

################################################################################

# make plots that show all exposure maps and inset anomally histograms
# saves 6-8 plots for each species to a folder -- to combine manually

group_exposure_plots <- function(number, this_path_name, this_species_layer, this_EFH_level, this_species_name){
  exposure_map_list = list()

  exposure_factors_list <- layer_names[number,6:13] # identify species in layer list 
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)]) # identify set of exposure factors for this species
  
  for(i in 1:length(these_exposure_factors$V1)){
    original_exposure_data <- create_overlap(this_path_name, this_species_layer, this_EFH_level, this_species_name, anomaly[[paste(these_exposure_factors$V1[i], "_anomaly", sep = "")]], these_exposure_factors$V1[i])
    original_exposure_data <- original_exposure_data |>
      dplyr::select(anomaly, geometry, type, anomaly_bins) |>
      rename(variable = type)
    exposure_map_list[[i]] <- original_exposure_data # add it to the list
  }
  
  # create data frame to put exposure score numbers in
  species_scores <- data.frame(
    exposure_factor = unlist(these_exposure_factors$V1)
  )
  
  for(i in 1:length(these_exposure_factors$V1)){
    # calculate exposure score for each exposure factor
    score <- calculate_exposure_score(original_exposure_data, this_path_name, this_species_layer, this_EFH_level, this_species_name, anomaly[[paste(these_exposure_factors$V1[i], "_anomaly", sep = "")]], these_exposure_factors$V1[i])
    species_scores[i,"score"] <- score
  }
  
  # group_exposure_maps <- do.call(rbind, exposure_map_list) 
  # group_exposure_plot <- do.call(rbind, exposure_plot_list) 
  
  exposure_maps <- function(map_data, score){
    map <- ggplot() + 
      geom_sf(data = map_data, aes(color = anomaly, geometry = geometry), size = 1, alpha = 0.8) +
      coord_sf(xlim = c(-170, -130), ylim = c(50, 62)) +
      # facet_wrap(~variable, ncol = 1) +
      geom_sf(data = GOA, size = 0.2, fill = "gray85") +
      geom_sf(data = canada, size = 0.2, fill = "gray95") +
      scale_x_continuous(n.breaks = 4, expand = c(0,0)) +
      scale_y_continuous(breaks = seq(52, 60, by = 4), expand = c(0,0)) +
      scale_color_gradientn(
        rescaler = function (...) {
          scales::rescale_mid(..., mid = 0)
        },
        colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
        values = scales::rescale(c(-2, -1.5, -0.5, 0.5, 1.5, 2)),
        name = "Anomaly") +
      labs(x = "Longitude",
           y = "Latitude",
           color = "Anomaly") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = c(0.95,0.85), legend.direction = "vertical", 
            legend.text = element_text(size = 6), legend.title = element_text(size = 8, margin = margin(b = 5)), 
            legend.key.size = unit(0.25, "cm"), legend.key.spacing = unit(0.05, "cm"),
            legend.frame = element_rect(color = "black", linewidth = 0.2),
            legend.background = element_rect(fill = "transparent"),
            legend.ticks = element_line(color = "black", linewidth = 0.25),
            strip.text = element_text(hjust = 0, size = 10),
            strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
            panel.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
            plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
      )
    
    inset <- ggplot(map_data) +
      geom_histogram(aes(x = anomaly, y = after_stat(count / sum(count)), fill = anomaly_bins), binwidth = 0.25, boundary = 0, linewidth = 0.25, colour="black", show.legend = FALSE) +
      scale_fill_manual(values = c("low" = "green", 
                                   "moderate" = "yellow", 
                                   "high" = "orange", 
                                   "very high" = "red")) +
      scale_x_continuous(expand = c(0,0)) +
      xlab("Anomaly") +
      ylab("Proportion") +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 5),
            rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
    
    ggdraw() +
      draw_plot(map) +
      draw_plot(inset, height = 0.15, x = 0.65, y = 0.15, vjust = 1)
    
    combined_plot <- map +
      inset_element(inset, align_to = "plot", left = 0.5, bottom = 0.17, right = 0.8, top = 0.57) 
    combined_plot$patches$layout$widths  <- 0.5
    combined_plot$patches$layout$heights <- 0.5
    
    return(combined_plot)
  }

#pollock <- exposure_maps(exposure_map_list[[1]], species_scores[1, "score"])
#ggsave(filename="pollock_bt.png", path = here("plots/"), plot = pollock, device = "png", width = 7, height = 4, bg = "transparent", dpi = 400)

exposure_plot_list = list()
for(i in 1:length(exposure_map_list)){
  data <- exposure_map_list[[i]]
  plot <- exposure_maps(data, species_scores[i, "score"])
  
  exposure_plot_list[[i]] <- plot
  
  ggsave(filename=paste("exposure", species_scores[i, "exposure_factor"], ".png", sep = ""), path = here(paste("plots/", this_species_name, sep = "")), plot = plot, device = "png", width = 7, height = 4, bg = "transparent", dpi = 400)
}

#plot_pollock_exposure <- patchwork::wrap_plots(exposure_plot_list, ncol = 2, axes = "collect", axis_titles = "collect")
#ggsave(filename="pollock_exposure_all.png", path = here("plots/"), plot = plot_pollock_exposure, device = "png", width = 7, height = 4, bg = "transparent", dpi = 400)
}

# run for pollock
group_exposure_plots(54, "gdb_path", "GOA_adult_walleyepollock_efh_level2_abundance_summer", "2", "Walleye pollock")
group_exposure_plots(44, "gdb_path", "GOA_adult_shortrakerrockfish_efh_level2_abundance_summer", "2", "Shortraker rockfish")
group_exposure_plots(17, "bts_sdm_path", "predictions_grenadiers.rda", "2", "Giant grenadier")
group_exposure_plots(34, "gdb_path", "GOA_adult_petralesole_efh_level2_abundance_summer", "2", "Petrale sole")
