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

exposure_factors <- c("SST", "BT", "salinity")

# ============================================================================
# SECTION 2: Calculate Anomaly
# ============================================================================
# Calculate standardized anomaly: (future mean - historical mean) / historical SD.
# ============================================================================

# create function that calculates an anomaly for all exposure factors in this script

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
BT_anomaly <- create_anomaly(temp_bottom_future, temp_bottom_hindcast)
salinity_anomaly <- create_anomaly(salt_surface_future, salt_surface_hindcast)

coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

create_anomaly_plot <- function(data, exposure_factor_name) {
  # Load coastline data for map visualization
  
  plot <- ggplot() + # plot climate anomaly
    geom_sf(data = data, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = coast, color = "black", linewidth = 0.3) +
    scale_color_gradientn(
      colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
      values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
      limits = c(-3, 3),
      name = "Anomaly"
    ) +
    labs(x = "Longitude",
         y = "Latitude",
         color = "standardized anomaly") +
    theme_bw() +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1.5)) +
    theme(legend.position = c(0.6, 0.15), legend.direction = "horizontal", 
          legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
          legend.key.size = unit(0.7, "cm"), legend.key.spacing = unit(0.12, "cm"),
          legend.frame = element_rect(color = "black", linewidth = 0.25),
          legend.background = element_rect(fill = "transparent"),
          legend.ticks = element_line(color = "black", linewidth = 0.25),
          strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          plot.background = element_rect(fill = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  print(plot)
  ggsave(paste(exposure_factor_name, "anomaly.png", sep="_"), path = here("plots/Anomaly plots"), device = "png", dpi = 300, width=8, height=5)

}

create_anomaly_plot(SST_anomaly, "SST")
create_anomaly_plot(BT_anomaly, "BT")
create_anomaly_plot(salinity_anomaly, "salinity")

# ============================================================================
# SECTION 3: Calculate exposure
# ============================================================================
# Map EFH maps and create overlap with ROMS outputs. Determine exposure and 
# create associated plots.
# ============================================================================

####### Load, Adjust, and Plot EFH Maps ######

#tell R where files are 
gdb_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/GOA_groundfish_2023.gdb"
scallop_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/EFH_2018_Scallop.gdb"
salmon_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/Salmon_2023.gdb"

#list all vector layers in the .gbd
layers <- vector_layers(gdb_path)
print(layers)

scallop_layers <- vector_layers(scallop_path)
print(scallop_layers)

salmon_layers <- vector_layers(salmon_path)
print(salmon_layers)

# get list of all the species layers we want to put through the function
layer_names <- read_excel(here("goa_efh_spp_lifestages.xlsx"), sheet = "groundfish",
                           col_names = c("species_name", "layer_name"))
species_layers<-unique(layer_names$layer_name)

# set colors for plotting EFH areas
EFH_cols <- c( 
  "2" = "#4B0082", #value <= 2
  "3" = "#006FA5", #value <= 3
  "4" = "#3CB371", #value <= 4
  "5" = "#FFFF00"  #value <= 5
)

create_EFH_layer <- function(species_layer){
  
  #load a specific layer
  filtered <- vect(gdb_path, layer = species_layer)  
  
  # adjust the EFH map 
  filtered <- project(filtered, "epsg:4326") # project EFH map onto same crs as ROMS 
  sf_filtered <- st_as_sf(filtered, crs = 4326) |> # make sf object
    filter(layer != "1") |> # remove non-EFH areas
    st_make_valid() # make geometry valid (not sure why)
  sf_use_s2(FALSE) # don't use spherical geometry (not sure why)
  sf_filtered <- sf_filtered |> st_shift_longitude() # Convert to 0-360° longitude to match ROMS data
  sf_filtered$layer <- as.character(sf_filtered$layer) # fix EFH layer class to text instead of numbers
  
  return(sf_filtered)
}

plot_EFH_layer <- function(species_layer){
  sf_filtered = create_EFH_layer(species_layer)
  plot <- ggplot() +
    geom_sf(data = sf_filtered, aes(color = layer, geometry = geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = coast, color = "black", linewidth = 0.3) +
    scale_color_manual(
      values = EFH_cols,
      labels = c(
        "2" = "95% EFH Area",
        "3" = "75% Principal EFH Area",
        "4" = "50% Core EFH Area",
        "5" = "25% EFH Hot Spots"
      ),
      # 2. Control the legend breaks and order
      breaks = c("2", "3", "4", "5")
    ) +
    labs(x = "Longitude",
         y = "Latitude",
         color = "EFH area") +
    theme_bw() +
    theme(legend.position = c(0.61, 0.2), legend.direction = "vertical", 
          legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
          legend.key.size = unit(0.5, "cm"), legend.key.spacing = unit(0.12, "cm"),
          legend.background = element_rect(fill = "transparent"),
          strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  ggsave(paste(species_layer, "EFH.png", sep="_"), device = "png", path = here("plots/EFH plots"), plot = plot, width=8, height=5, dpi = 300)
}

# testing
plot_EFH_layer("GOA_adult_alaskaplaice_efh_level2_abundance_summer")

# for loop to create EFH plots for all species (30) 
for (i in 1:length(species_layers)) {
  plot_EFH_layer(gdb_path, species_layers[i])
}

####### Calculate exposure ######

# function to calculate exposure for EFH species
create_overlap <- function(species_layer, species_name, anomaly_data, exposure_factor_name){
  sf_filtered <- create_EFH_layer(species_layer)
  
  anomaly_data = st_transform(anomaly_data, crs = st_crs(sf_filtered)) # match crs of ROMS points and EFH polygons ??
  
  # find intersects between points (ROMS) and polygons (EFH)
  anomaly_data$EFH <- apply(st_intersects(sf_filtered, anomaly_data, sparse = FALSE), 2, 
                             function(col) {sf_filtered[which(col), ]$layer}) # not exactly sure how this code works, but it seems to have done what I want it to
  
  # remove every point that lies outside of the species distribution
  exposure <- anomaly_data |> filter(EFH == "4" | EFH == "5") # filter to core habitat
  
  #dir.create(here(paste("plots/Exposure plots/", species_name, "/", exposure_factor_name, sep ="")), recursive = TRUE)
  return(exposure)
}

assign_exposure_levels <- function(species_layer, species_name, anomaly_data, exposure_factor_name){
  exposure <- create_overlap(species_layer, species_name, anomaly_data, exposure_factor_name)
  # assign categories from low - very high to the anomaly values
  exposure_plot <- exposure |> 
    mutate(
      exposure_score = ifelse(anomaly >= -0.5 & anomaly <= 0.5, "low", ifelse((anomaly < -0.5 & anomaly >= -1.5) | (anomaly > 0.5 & anomaly <= 1.5), "moderate", ifelse((anomaly < -1.5 & anomaly >= -2) | (anomaly > 1.5 & anomaly <= 2), "high", "very_high")))
    ) |>
    st_drop_geometry()
  
  # set levels for exposure scores from low - very high
  exposure_plot$exposure_score <- factor(exposure_plot$exposure_score, levels = c("low", "moderate", "high", "very_high"))
  exposure_plot$exposure_score = ordered(exposure_plot$exposure_score,
                                         levels = c("low",
                                                    "moderate",
                                                    "high",
                                                    "very_high"))
  
  # calculate counts and proportion in each scoring category 
  exposure_plot <- exposure_plot |>
    group_by(exposure_score, .drop = FALSE) |>
    summarize(count = n()) |>
    ungroup() |>
    complete(exposure_score, fill = list(count = 0)) |>
    mutate(
      total = sum(count),
      prop = count / total
    )
  return(exposure_plot)
}
  
calculate_exposure_score <- function(species_layer, species_name, anomaly_data, exposure_factor_name){
  exposure_plot <- assign_exposure_levels(species_layer, species_name, anomaly_data, exposure_factor_name)
  # calculate weighted mean
  exp_fact_mean <- exposure_plot |>
    select(!c(prop,total)) |>
    pivot_wider(names_from = "exposure_score",
                values_from = "count")|>
    mutate(
      weighted_mean = ((low*1)+(moderate*2)+(high*3)+(very_high*4))/(low+moderate+high+very_high)
    )
  exp_fact_mean$weighted_mean <- round(exp_fact_mean$weighted_mean, digits = 2)
  
  return(exp_fact_mean$weighted_mean)
}

create_exposure_plots <- function(species_layer, species_name, anomaly_data, exposure_factor_name) {
  original_exposure_data <- create_overlap(species_layer, species_name, anomaly_data, exposure_factor_name)
  exposure_plot <- assign_exposure_levels(species_layer, species_name, anomaly_data, exposure_factor_name)
  # plot exposure map
  plot2 <- ggplot() + 
    geom_sf(data = original_exposure_data, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = coast, color = "black", linewidth = 0.3) +
    scale_color_gradientn(
      colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
      values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
      limits = c(-3, 3)
    ) +
    labs(x = "Longitude",
         y = "Latitude",
         color = "Anomaly") +
    theme_bw() +
    theme(
      legend.position = c(0.6, 0.15), legend.direction = "horizontal", 
      legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
      legend.key.size = unit(0.7, "cm"), legend.key.spacing = unit(0.12, "cm"),
      legend.frame = element_rect(color = "black", linewidth = 0.25),
      legend.background = element_rect(fill = "transparent"),
      legend.ticks = element_line(color = "black", linewidth = 0.25),
      strip.text = element_text(hjust = 0, size = 10),
      strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
      plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
    )
  ggsave(paste(species_layer, "exposure_map.png", sep="_"), path = here(paste("plots/Exposure plots/", species_name, "/", exposure_factor_name, sep ="")), plot = plot2, device = "png", width = 8, height = 5, dpi = 300)
  
  # make histogram of anomalies
  plot3 <- ggplot(original_exposure_data) +
    geom_histogram(aes(x = anomaly, y = after_stat(count / sum(count)), fill = anomaly_bins), binwidth = 0.25, boundary = 0, linewidth = 0.25, colour="black", show.legend = FALSE) +
    scale_fill_manual(values = c("low" = "green", 
                                 "moderate" = "yellow", 
                                 "high" = "orange", 
                                 "very high" = "red")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(
      limits = c(-3,3),
      breaks = seq(-3, 3, by = 1), 
      expand = c(0,0)
    ) +
    xlab("Anomaly") +
    ylab("Percent") +
    theme_bw() +
    theme(rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  ggsave(paste(species_layer, exposure_factor_name, "anomaly_histogram.png", sep="_"), path = here(paste("plots/Exposure plots/", species_name, "/", exposure_factor_name, sep ="")), plot = plot3, device = "png")
  
  exp_fact_mean <- calculate_exposure_score(species_layer, species_name, anomaly_data, exposure_factor_name)
  
  # plot distribution of exposure scores
  plot4 <- ggplot(exposure_plot) +
    geom_col(mapping = aes(x = exposure_score, y = prop, fill = exposure_score), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
    scale_x_discrete(labels = c("L", "M", "H", "V")) +
    scale_y_continuous(labels = scales::percent) +
    ylab("Percent") +
    xlab("Exposure Score") +
    scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
    annotate("text", x = I(0.8), y = I(0.8), label = paste("Exposure = ", exp_fact_mean, sep = "")) + # add exposure score in top right corner
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  ggsave(paste(species_layer, exposure_factor_name, "exposure_scores.png", sep="_"), path = here(paste("plots/Exposure plots/", species_name, "/", exposure_factor_name, sep ="")), plot = plot4, device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)
  
}
  
species_name<-unique(layer_names$species_name)

# tests - ak plaice
assign_exposure_levels("GOA_adult_alaskaplaice_efh_level2_abundance_summer", "Alaska plaice", SST_anomaly, "SST")
calculate_exposure_score("GOA_adult_alaskaplaice_efh_level2_abundance_summer", "Alaska plaice", SST_anomaly, "SST")
create_exposure_plots("GOA_adult_alaskaplaice_efh_level2_abundance_summer", "Alaska plaice", SST_anomaly, "SST")

###### Create exposure plots - for loops by exposure factor #######

# for loop to create all species (30) exposure plots for SST
for (i in 1:length(species_layers)) {
  create_exposure_plots(species_layers[i], species_name[i], SST_anomaly, "SST") 
}
# for loop to create all species (30) exposure plots for BT
for (i in 1:length(species_layers)) {
  create_exposure_plots(species_layers[i], species_name[i], BT_anomaly, "BT")
}
# for loop to create all species (30) exposure plots for salinity
for (i in 1:length(species_layers)) {
  create_exposure_plots(species_layers[i], species_name[i], salinity_anomaly, "salinity")
}

###### Calculate exposure scores - for loops by exposure factor #######

# for loop to create all species (30) exposure plots for SST
for (i in 1:length(species_layers)) {
  exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], SST_anomaly, "SST")
  layer_names[i,"SST"] <- exposure_score
}
# for loop to create all species (30) exposure plots for BT
for (i in 1:length(species_layers)) {
  exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], BT_anomaly, "BT")
  layer_names[i,"BT"] <- exposure_score
}
# for loop to create all species (30) exposure plots for salinity
for (i in 1:length(species_layers)) {
  exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], salinity_anomaly, "salinity")
  layer_names[i,"salinity"] <- exposure_score
}

# original code (for just SST):
#   exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], SST_anomaly, "SST")
#   layer_names[i,"SST] <- exposure_score

# # make for loop to calculate exposure scores for all listed exposure factors then add them to the table layer_names
# for(j in 1:length(exposure_factors)) {
# for (i in 1:length(species_layers)) {
#   exposure_score <- calculate_exposure_score(species_layers[i], species_name[i], SST_anomaly, exposure_factors[j])
#   layer_names[i,exposure_factors[j]] <- exposure_score
# }
# }

