# ============================================================================
# Draft exposure script
# ============================================================================
# Purpose: This script calculate exposure for a single species (AK plaice or Walleye pollock) and
# single exposure factor (sea surface temperature). It produces 3 plots:
#   1. Map of SST anomalies across the AK plaice or Walleye pollock distribution in the GOA.
#   2. Histogram of anomalies across AK plaice or Walleye pollock distribution.
#   3. Distribution of exposure scores and ultimate exposure score.

# Clear work space and free up memory
rm(list = ls())
gc()

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, here, terra)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
here::i_am("exposure_example.R")

# Load custom analysis functions by Albi
source(here("functions.R"))

# ============================================================================
# SECTION 1: Load and Filter Data - Historical and Future
# ============================================================================
# The Arrow package allows "lazy evaluation" - you can filter the data before
# loading it into RAM, which is much faster than loading everything first.
# ============================================================================

# Sea Surface Temperature
SST <- open_dataset(here("data/processed/all_scenarios_bias_corrected_temp_surface.parquet")) # load projections
SST_future <- SST |> filter(run == "ssp585") |> # future projections
  filter(date > as.Date("2030-01-01")) |>  # Restrict dates to later than 2030
  filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
  collect()  # NOW load the filtered data into RAM
SST_hindcast <- SST |> filter(run == "hindcast") |> # hindcast 
  filter(date > as.Date("1991-01-01")) |>  # Restrict dates to later than 1991
  filter(date < as.Date("2020-12-31")) |>  # Restrict to 2020
  collect()  # NOW load the filtered data into RAM

# Bottom Temperature
BT <- open_dataset(here("data/processed/all_scenarios_bias_corrected_temp_bottom.parquet")) # load projections

BT_future <- BT |> filter(run == "ssp585") |> # future projections
  filter(date > as.Date("2030-01-01")) |>  # Restrict dates to later than 2030
  filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
  collect()  # NOW load the filtered data into RAM
BT_hindcast <- BT |> filter(run == "hindcast") |> # hindcast 
  filter(date > as.Date("1991-01-01")) |>  # Restrict dates to later than 1991
  filter(date < as.Date("2020-12-31")) |>  # Restrict to 2020
  collect()  # NOW load the filtered data into RAM

# ============================================================================
# SECTION 2: Calculate Sea Surface and Bottom Temperature Anomalies
# ============================================================================
# Calculate standardized anomaly: (future mean - historical mean) / historical SD 
# for SST and BT.
# ============================================================================

# sea surface tempeature
data_anomaly <- SST_future |> # calculate climate anomaly (future mean - historical mean / historical standard deviation)
  filter(month == "7" | month == "8" | month == "9") |> # can filter by winter or summer months by changing month numbers
  summarize(average_future = mean(value_dc), .by = c(cell_id, lon_rho, lat_rho)) |>
  left_join(
    SST_hindcast |> 
      filter(month == "7" | month == "8" | month == "9") |> # filter in hindcast as well
      summarize(average_hist = mean(value), sd_hist = sd(value), .by = c(cell_id, lon_rho, lat_rho)), # for hindcast, value=value_dc
    by = join_by(lon_rho, lat_rho)
  ) |>
  mutate(anomaly = (average_future-average_hist)/sd_hist) |> # calculate anomaly 
  select(!c(average_future, average_hist, sd_hist))

# Load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

temp_anoms <- data_anomaly |> # set scoring categories
  st_as_sf(coords = c("lon_rho", "lat_rho"))
st_crs(temp_anoms)= 4326
temp_anoms <- temp_anoms |>
  mutate( # set scoring categories
    anomaly_bins = case_when(anomaly >= -5 & anomaly < -2 ~"very high",
                             anomaly >= -2 & anomaly < -1.5 ~"high",
                             anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                             anomaly >= -0.5 & anomaly < 0.5~"low",
                             anomaly >= 0.5 & anomaly < 1.5~"moderate",
                             anomaly >= 1.5 & anomaly < 2~"high",
                             anomaly >=2 & anomaly <= 5~"very high")
  )

plot <- ggplot() + # plot climate anomaly
  geom_sf(data = temp_anoms, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_gradientn(
    colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
    values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
    limits = c(-3, 3),
    name = "Anomaly"
  ) +
  labs(title = "SST anomaly",
       x = "Longitude",
       y = "Latitude",
       color = "standardized anomaly") +
  theme_bw() +
  theme(legend.position = "bottom")
print(plot)
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="SSTanomaly.png", device = "png")

# bottom temperature (same process as above)
data_anomaly_bt <- BT_future |> # calculate climate anomaly (future mean - historical mean / historical standard deviation)
  filter(month == "7" | month == "8" | month == "9") |> # can filter by winter or summer months by changing month numbers
  summarize(average_future = mean(value_dc), .by = c(cell_id, lon_rho, lat_rho)) |>
  left_join(
    BT_hindcast |> 
      filter(month == "7" | month == "8" | month == "9") |> # filter in hindcast as well
      summarize(average_hist = mean(value), sd_hist = sd(value), .by = c(cell_id, lon_rho, lat_rho)), # for hindcast, value=value_dc
    by = join_by(lon_rho, lat_rho)
  ) |>
  mutate(anomaly = (average_future-average_hist)/sd_hist) |> # calculate anomaly 
  select(!c(average_future, average_hist, sd_hist))

# Load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

temp_anoms_bt <- data_anomaly_bt |> # set scoring categories
  st_as_sf(coords = c("lon_rho", "lat_rho"))
st_crs(temp_anoms_bt)= 4326
temp_anoms_bt <- temp_anoms_bt |>
  mutate( # set scoring categories
    anomaly_bins = case_when(anomaly >= -5 & anomaly < -2 ~"very high",
                             anomaly >= -2 & anomaly < -1.5 ~"high",
                             anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                             anomaly >= -0.5 & anomaly < 0.5~"low",
                             anomaly >= 0.5 & anomaly < 1.5~"moderate",
                             anomaly >= 1.5 & anomaly < 2~"high",
                             anomaly >=2 & anomaly <= 5~"very high")
  )

plot <- ggplot() + # plot climate anomaly
  geom_sf(data = temp_anoms_bt, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_gradientn(
    colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
    values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
    limits = c(-3, 3),
    name = "Anomaly"
  ) +
  labs(title = "BT anomaly",
       x = "Longitude",
       y = "Latitude",
       color = "standardized anomaly") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
print(plot)
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="BTanomaly.png", device = "png")

# ============================================================================
# SECTION 3: Plot species distribution 
# ============================================================================
# Map EFH maps and creating overlap with ROMS outputs.
# ============================================================================

#tell R where files are 
gdb_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/GOA_groundfish_2023.gdb"

#list all vector layers in the .gbd
layers <- vector_layers(gdb_path)
print(layers)

#load a specific layer
adult_alaskaplaice <- vect(gdb_path, layer = "GOA_adult_alaskaplaice_efh_level2_abundance_summer")
adult_walleyepollock <- vect(gdb_path, layer = "GOA_adult_walleyepollock_efh_level2_abundance_summer")

# set colors for plotting
EFH_cols <- c( 
  "2" = "#4B0082", #value <= 2
  "3" = "#006FA5", #value <= 3
  "4" = "#3CB371", #value <= 4
  "5" = "#FFFF00"  #value <= 5
)

# adjust the EFH map - AK plaice
adult_alaskaplaice <- project(adult_alaskaplaice, "epsg:4326") # project EFH map onto same crs as ROMS 
sf_akplaice <- st_as_sf(adult_alaskaplaice, crs = 4326) |> # make sf object
  filter(layer != "1") |> # remove non-EFH areas
  st_make_valid() # make geometry valid (not sure why)
sf_use_s2(FALSE) # don't use spherical geometry (not sure why)
sf_akplaice <- sf_akplaice |> st_shift_longitude() # Convert to 0-360° longitude to match ROMS data
sf_akplaice$layer <- as.character(sf_akplaice$layer) # fix EFH layer class to text instead of numbers

# adjust the EFH map - walleye pollock
adult_walleyepollock <- project(adult_walleyepollock, "epsg:4326") # project EFH map onto same crs as ROMS 
sf_walleyepollock <- st_as_sf(adult_walleyepollock, crs = 4326) |> # make sf object
  filter(layer != "1") |> # remove non-EFH areas
  st_make_valid() # make geometry valid (not sure why)
sf_use_s2(FALSE) # don't use spherical geometry (not sure why)
sf_walleyepollock <- sf_walleyepollock |> st_shift_longitude() # Convert to 0-360° longitude to match ROMS data
sf_walleyepollock$layer <- as.character(sf_walleyepollock$layer) # fix EFH layer class to text instead of numbers

# plot the EFH map
ggplot() +
  geom_sf(data = sf_walleyepollock, aes(color = layer, geometry = geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_manual(
    values = EFH_cols,
    labels = c(
      "2" = "95% EFH Area (All shaded areas)",
      "3" = "75% Principal EFH Area",
      "4" = "50% Core EFH Area",
      "5" = "25% EFH Hot Spots"
    ),
    # 2. Control the legend breaks and order
    breaks = c("2", "3", "4", "5")
  ) +
  labs(title = "EFH Map",
       x = "Longitude",
       y = "Latitude",
       color = "EFH area") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
ggsave(filename="walleyepollock_EFH.png", device = "png")

temp_anoms_bt = st_transform(temp_anoms_bt, crs = st_crs(sf_walleyepollock)) # match crs of ROMS points and EFH polygons ??

# find intersects between points (ROMS) and polygons (EFH)
temp_anoms_bt$EFH <- apply(st_intersects(sf_walleyepollock, temp_anoms_bt, sparse = FALSE), 2, 
                        function(col) {sf_walleyepollock[which(col), ]$layer}) # not exactly sure how this code works, but it seems to have done what I want it to

# remove every point that lies outside of the species distribution
exposure <- temp_anoms_bt |> filter(EFH == "2" | EFH == "3" | EFH == "4" | EFH == "5")

plot <- ggplot() + # plot exposure map
  geom_sf(data = exposure, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_gradientn(
    colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
    values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
    limits = c(-3, 3)
  ) +
  labs(title = "Exposure Map",
       x = "Longitude",
       y = "Latitude",
       color = "Anomaly") +
  theme_bw() +
  theme(
    strip.text = element_text(hjust = 0, size = 10),
    strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
    rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="walleyepollock_bt_exposure map.png", device = "png")

# ============================================================================
# SECTION 4: Calculate exposure scores  
# ============================================================================
# Apply logic model and calculate exposure score. Create histogram summaries of 
# anomalies and exposure scores.
# ============================================================================

# make histogram of anomalies
ggplot(exposure) +
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
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="walleyepollock_bt_anomaly_histogram.png", device = "png")

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

# calculate weighted mean
exp_fact_mean <- exposure_plot |>
  select(!c(prop,total)) |>
  pivot_wider(names_from = "exposure_score",
               values_from = "count")|>
  mutate(
    weighted_mean = ((low*1)+(moderate*2)+(high*3)+(very_high*4))/(low+moderate+high+very_high)
  )
exp_fact_mean$weighted_mean <- round(exp_fact_mean$weighted_mean, digits = 2)

# plot distribution of exposure scores
ggplot(exposure_plot) +
  geom_col(mapping = aes(x = exposure_score, y = prop, fill = exposure_score), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
  scale_x_discrete(labels = c("L", "M", "H", "V")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent") +
  xlab("Exposure Score") +
  scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
  annotate("text", x = I(0.8), y = I(0.8), label = paste("Exposure = ", exp_fact_mean$weighted_mean, sep = "")) + # add exposure score in top right corner
  theme_bw() +
  theme(strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="walleyepollock_bt_exposure_scores.png", plot = get_last_plot(), device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)

