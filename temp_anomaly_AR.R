# ============================================================================
# ROMS Temperature Anomaly Script
# ============================================================================
# Purpose: This script attempts to create a standardized climate anomalies for ROMS
# exposure variables in the Gulf of Alaska.
#
# What this script does:
#   1. Opens a processed ROMS data set (stored as Parquet) corresponding to 
#      future and hindcast projections
#   2. Filters to variable (e.g., temperature), depth, and time period 
#   3. Calculates temperature anomaly 
#   4. Creates a spatial map of the temperature anomaly
#
# IMPORTANT NOTES:
# - The Parquet files are pre-filtered to the Gulf of Alaska region and 
#   depths ≤1000m to minimize file size
# - If you need broader spatial coverage or deeper depths, you must reprocess
#   the original netCDF files using the process_annual_file() function
# - Projection runs are NOT bias-corrected yet
# - This script assumes you have already processed raw ROMS outputs into
#   Parquet format using the main processing workflow
# ============================================================================

# Clear workspace and free up memory
rm(list = ls())
gc()

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data

# Load custom analysis functions by Albi
source("functions.R")

# ============================================================================
# SECTION 1: Load and Filter Data - Historical and Future
# ============================================================================
# The Arrow package allows "lazy evaluation" - you can filter the data before
# loading it into RAM, which is much faster than loading everything first.
# ============================================================================

# # Load projections (not bias corrected)
# 
# # Specify which model run to analyze
# # Options  include: "hindcast", "historical", "ssp126", "ssp245", "ssp585"
# # WARNING: Projection runs (ssp*) are not yet bias-corrected
# this_run <- "ssp585"
# 
# # Construct the file path to the processed Parquet file
# this_file <- paste0("data/processed/", this_run, "_annual_data.parquet")
# 
# # Open the dataset WITHOUT loading it into memory
# # This creates a connection that allows you to query the data structure
# all_data <- open_dataset(this_file)  # Uses Apache Arrow for efficient access
# 
# # Filter to specific data subset before loading
# # This is much more memory-efficient than loading everything first
# # You don't have to apply filters but it will occupy less RAM if you do
# # If you want to read many data stream at once, you should definitely filter
# data_subset_future <- all_data %>%
#   filter(layer == "surface") %>%  # Select only surface layer (other options: "bottom")
#   filter(variable == "temp") %>%  # Choose variable of interest ("temp", "salt", "PhS", "PhL", etc.)
#   filter(date > as.Date("2030-01-01")) %>%  # Restrict dates to later than 2030
#   filter(date < as.Date("2059-12-31")) %>%  # Restrict to earlier than 2039
#   collect()  # NOW load the filtered data into RAM
# 
# # Load hindcasts (not bias corrected)
# 
# # Specify which model run to analyze
# # Options  include: "hindcast", "historical", "ssp126", "ssp245", "ssp585"
# # WARNING: Projection runs (ssp*) are not yet bias-corrected
# this_run <- "hindcast"
# 
# # Construct the file path to the processed Parquet file
# this_file <- paste0("data/processed/", this_run, "_annual_data.parquet")
# 
# # Open the dataset WITHOUT loading it into memory
# # This creates a connection that allows you to query the data structure
# all_data <- open_dataset(this_file)  # Uses Apache Arrow for efficient access
# 
# # Filter to specific data subset before loading
# # This is much more memory-efficient than loading everything first
# # You don't have to apply filters but it will occupy less RAM if you do
# # If you want to read many data stream at once, you should definitely filter
# data_subset_hist <- all_data %>%
#   filter(layer == "surface") %>%  # Select only surface layer (other options: "bottom")
#   filter(variable == "temp") %>%  # Choose variable of interest ("temp", "salt", "PhS", "PhL", etc.)
#   filter(date > as.Date("1991-01-01")) %>%  # Restrict dates to later than 1991
#   filter(date < as.Date("2020-12-31")) %>%  # Restrict to 2020
#   collect()  # NOW load the filtered data into RAM

# Bias corrected SST

SST <- open_dataset(here("data/processed/all_scenarios_bias_corrected_temp_surface.parquet")) # load projections
SST_future <- SST |> filter(run == "ssp585") |> # future projections
  filter(date > as.Date("2030-01-01")) |>  # Restrict dates to later than 2030
  filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
  collect()  # NOW load the filtered data into RAM
SST_hindcast <- SST |> filter(run == "hindcast") |> # hindcast 
  filter(date > as.Date("1991-01-01")) |>  # Restrict dates to later than 1991
  filter(date < as.Date("2020-12-31")) |>  # Restrict to 2020
  collect()  # NOW load the filtered data into RAM

# ============================================================================
# SECTION 2: Calculate Temperature Anomaly
# ============================================================================
# Calculate standardized anomaly: (future mean - historical mean) / historical SD 
# for SST.
# ============================================================================

# Notes AR:
# You are using the correct bias corrected projections, BUT - you should use "value_dc" for your calculations, as that is the delta-corrected value. For the hindcast, value=value_dc
# You take annual averages for temp. However, focusing on either summer or winter will give you more contrast in your anomalies. In the last step here you divide by the SD of the hindcast. There is a lot of variability over one year so your SD in the denominator is going to swamp the anomalies at the numerators if you took SD over the whole year
# You use a pretty long period as "future". I understand that this is meant to absorb short-term variability and I remember Al recommending using multi-decadal windows. However, hindcast and future period are pretty close to eachother in time and this is also going to dampen your anomalies
# Your future period is 2030-2059. ROMS really start diverging past 2060 as far as I can tell, so you should not expect extreme increases in temp before then
# I suggest you keep cell_id in the grouping. lon and lat get rounded sometimes you end up lumping a few cells together, which adds to the SD
# note that this is surface temperature so it may or may not be appropriate for groundfish

# see below for what happens if you use value_dc, summer (July-Sept) temperatures for both proj and hind (and keep cell_id)

# but first, check out these couple plots on the raw temp values, they should give you some insights
SST_checks <- SST |> filter(run %in% c("hindcast", "ssp126", "ssp585")) |> # future projections
  filter(date > as.Date("1991-01-01")) |>  # Restrict dates to later than 2030
  filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
  collect()  # NOW load the filtered data into RAM

# first, check out the difference between value and value_dc 
# I am using one month for ease of visualization - picking january because I know that the bias in GFLD is largest in the winter. You can fiddle with other months to learn something about the bias
SST_checks %>%
  filter(run != "ssp126", month == 1) %>%
  group_by(year) %>%
  summarise(mean_value = mean(value),
            mean_value_dc = mean(value_dc)) %>%
  pivot_longer(-year) %>%
  ggplot(aes(x = year, y = value, color = name))+geom_line()
  
# now look at how noisy the time series is if you use all months to get your annual averages - all that noise goes into your SD and swamps your anomalies
SST_checks %>%
  filter(run != "ssp126") %>%
  # filter(month %in% c(7:9)) %>% # toggle this line on and off to see how your time series changes as you filter by month
  group_by(year, month, run) %>%
  summarise(mean_value_dc = mean(value_dc)) %>%
  ungroup() %>%
  mutate(date = as.Date(sprintf("%d-%02d-15", year, month))) %>%
  dplyr::select(date, run, mean_value_dc) %>%
  ggplot(aes(x = date, y = mean_value_dc, color = run))+geom_line()+geom_point()

# now look how close ssp126 and ssp585 are in 2030-2059
SST_checks %>%
  filter(month %in% c(7:9)) %>% 
  group_by(year, run) %>%
  summarise(mean_value_dc = mean(value_dc)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = mean_value_dc, color = run))+geom_line()+geom_point()

# and here is your code using value_dc for the future, summer temps, and retaining cell_id
data_anomaly <- SST_future |> # calculate climate anomaly (future mean - historical mean / historical standard deviation)
  filter(month %in% 7:9) |> # can filter by winter or summer months by changing month numbers
  summarize(average_future = mean(value_dc), .by = c(cell_id, lon_rho, lat_rho)) |>
  left_join(
    SST_hindcast |> 
      filter(month %in% 7:9) |> # needs to be consistent between hidcast and future
      summarize(average_hist = mean(value), sd_hist = sd(value), .by = c(cell_id, lon_rho, lat_rho)), # for hindcast, value=value_dc
    by = join_by(cell_id, lon_rho, lat_rho)
    ) |>
  mutate(anomaly = (average_future-average_hist)/sd_hist) # calculate anomaly

# Load coastline data for map visualization
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

plot_data <- data_anomaly |> # set scoring categories
  st_as_sf(coords = c("lon_rho", "lat_rho"), crs = 4326) |>
  mutate(
    anomaly_bins = case_when(anomaly >= -5 & anomaly < -2 ~"very high",
                             anomaly >= -2 & anomaly < -1.5 ~"high",
                             anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                             anomaly >= -0.5 & anomaly < 0.5~"low",
                             anomaly >= 0.5 & anomaly < 1.5~"moderate",
                             anomaly >= 1.5 & anomaly < 2~"high",
                             anomaly >=2 & anomaly <= 5~"very high")
  )

plot <- ggplot() + # plot climate anomaly
  geom_sf(data = plot_data, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_gradientn(
    colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
    values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
    limits = c(-3, 3),
    name = "Anomaly"
  ) +
  labs(title = "Exposure Map",
       x = "Longitude",
       y = "Latitude",
       color = "standardized anomaly") +
  theme_bw() +
  theme(legend.position = "bottom")
print(plot)
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="SSTanomaly.png", device = "png")
  
# ============================================================================
# SECTION 3: Calculate exposure scores  
# ============================================================================
# Apply logic model and calculate exposure score. Create histogram summaries of 
# anomalies and exposure scores.
# ============================================================================
  
# make histogram of anomalies
ggplot(plot_data) +
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
  xlab("SST anomaly") +
  ylab("Percent") +
  theme_bw() +
  theme(rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="anomaly_histogram.png", device = "png")
  
# assign categories from low - very high to the anomaly values
exposure_plot <- data_anomaly |> 
  mutate(
    exposure_score = ifelse(anomaly >= -0.5 & anomaly <= 0.5, "low", ifelse((anomaly < -0.5 & anomaly >= -1.5) | (anomaly > 0.5 & anomaly <= 1.5), "moderate", ifelse((anomaly < -1.5 & anomaly >= -2) | (anomaly > 1.5 & anomaly <= 2), "high", "very high")))
  )

# set levels for exposure scores from low - very high
exposure_plot$exposure_score <- factor(exposure_plot$exposure_score, levels = c("low", "moderate", "high", "very high"))
exposure_plot$exposure_score = ordered(exposure_plot$exposure_score,
                                                 levels = c("low",
                                                            "moderate",
                                                            "high",
                                                           "very high"))

# calculate counts and proportion in each scoring category 
exposure_plot <- exposure_plot |>
  group_by(exposure_score, .drop =FALSE) |>
  summarize(count = n()) |>
  ungroup() |>
  complete(exposure_score,fill = list(count = 0)) |>
  mutate(
    total = sum(count),
    prop = count / total
  )

# plot distribution of exposure scores
ggplot(exposure_plot) +
  geom_col(mapping = aes(x = exposure_score, y = prop, fill = exposure_score), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
  scale_x_discrete(labels = c("L", "M", "H", "V")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent") +
  xlab("Exposure Score") +
  scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
setwd("~/Documents/OSU/GOA CVA/Exposure/CVA/plots/")
ggsave(filename="exposure_scores.png", plot = get_last_plot(), device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)

