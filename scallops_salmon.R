# Scallops and Salmon (presence/absence)

# Load required packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidync, ncdf4, lubridate, here, sf, rnaturalearth, arrow, dplyr, here, terra)
pacman::p_load_gh("ropensci/rnaturalearthhires")  # High-resolution coastline data
here::i_am("scallops_salmon.R")

scallop_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/EFH_2018_Scallop.gdb"
salmon_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/Salmon_2023.gdb"

scallop_layers <- vector_layers(scallop_path)
print(scallop_layers)

salmon_layers <- vector_layers(salmon_path)
print(salmon_layers)

# coastline for plots
coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%  # Crop to GOA region
  st_shift_longitude()  # Convert to 0-360° longitude to match ROMS data

######### SCALLOP ##########

#load scallop layer
scallop_EFH <- vect(scallop_path, layer = "_Weathervane_scallop_adult_EFH_Level1")  

# adjust the EFH map 
scallop_EFH <- project(scallop_EFH, "epsg:4326") # project EFH map onto same crs as ROMS 
sf_scallop_EFH <- st_as_sf(scallop_EFH, crs = 4326) |> # make sf object
  st_make_valid() # make geometry valid (not sure why)
sf_use_s2(FALSE) # don't use spherical geometry (not sure why)
sf_scallop_EFH <- sf_scallop_EFH |> st_shift_longitude() # Convert to 0-360° longitude to match ROMS data
# sf_scallop_EFH$layer <- as.character(sf_scallop_EFH$layer) # fix EFH layer class to text instead of numbers


ggplot() +
  geom_sf(data = sf_scallop_EFH, aes(geometry = geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(legend.position = c(0.61, 0.2), legend.direction = "vertical", 
        legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
        legend.key.size = unit(0.5, "cm"), legend.key.spacing = unit(0.12, "cm"),
        legend.background = element_rect(fill = "transparent"),
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))

ggsave("scallop_EFH.png", device = "png", path = here("plots/EFH plots"), plot = get_last_plot(), width=8, height=5, dpi = 300)

# Sea Surface Temperature
SST <- open_dataset(here("data/processed/all_scenarios_bias_corrected_temp_surface.parquet")) # load projections
SST_future <- SST |> filter(run == "ssp585") |> # future projections
  filter(date > as.Date("2030-01-01")) |>  # Restrict dates to later than 2030
  filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
  collect()  # NOW load the filtered data into RAM
SST_hindcast <- SST |> filter(run == "hindcast") |> # hindcast 
  filter(date > as.Date("1993-01-01")) |>  # Restrict dates to later than 1991
  filter(date < as.Date("2019-12-31")) |>  # Restrict to 2020
  collect()  # NOW load the filtered data into RAM
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

temp_anoms = st_transform(temp_anoms, crs = st_crs(sf_scallop_EFH)) # match crs of ROMS points and EFH polygons ??

# find intersects between points (ROMS) and polygons (EFH)
temp_anoms$EFH <- apply(st_intersects(sf_scallop_EFH, temp_anoms, sparse = FALSE), 2, 
                          function(col) {sf_scallop_EFH[which(col), ]$geometry}) # not exactly sure how this code works, but it seems to have done what I want it to

# remove every point that lies outside of the species distribution
exposure <- temp_anoms |> mutate(na=map_lgl(.x = EFH, .f = is_empty)) # identify rows not not included in the scallop polygon
exposure <- exposure |> filter(na == FALSE)

ggplot() + # plot exposure map for scallops
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
ggsave(filename="scallop_sst_exposure_map.png", device = "png", path = here("plots/Exposure plots/Weathervane scallop"))

######### SALMON ##########

sockeye_EFH <- vect(salmon_path, layer = "Salmon_mature_Sockeye_efh_level1_distribution")

# adjust the EFH map 
sockeye_EFH <- project(sockeye_EFH, "epsg:4326") # project EFH map onto same crs as ROMS 
sf_sockeye_EFH <- st_as_sf(sockeye_EFH, crs = 4326) |> # make sf object
  st_make_valid() # make geometry valid (not sure why)
sf_use_s2(FALSE) # don't use spherical geometry (not sure why)
sf_sockeye_EFH <- sf_sockeye_EFH |> st_shift_longitude() # Convert to 0-360° longitude to match ROMS data

# plot sockeye salmon EFH distribution
ggplot() +
  geom_sf(data = sf_sockeye_EFH, aes(geometry = geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  labs(x = "Longitude",
       y = "Latitude",) +
  theme_bw() +
  theme(legend.position = c(0.61, 0.2), legend.direction = "vertical", 
        legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
        legend.key.size = unit(0.5, "cm"), legend.key.spacing = unit(0.12, "cm"),
        legend.background = element_rect(fill = "transparent"),
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))

ggsave("sockeye_EFH.png", device = "png", path = here("plots/EFH plots"), plot = get_last_plot(), width=8, height=5, dpi = 300)


