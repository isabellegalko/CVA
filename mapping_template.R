# ============================================================================
# EFH Mapping Script
# ============================================================================
# Purpose: Mapping EFH maps and creating overlap with ROMS outputs
# Author: Mason Smith, Isabelle Galko
# Date created: 8 April 2026
# Date modified: 14 April 2026

#load packages
library(terra) #vector and raster files
library(ggplot2) #plotting
library(sf)

#tell R where files are 
gdb_path <- "/Users/isabellegalko/Documents/OSU/GOA CVA/Exposure/CVA/data/GOA_groundfish_2023.gdb"

#list all vector layers in the .gbd
layers <- vector_layers(gdb_path)
print(layers)

#load a specific layer
adult_alaskaplaice <- vect(gdb_path, layer = "GOA_adult_alaskaplaice_efh_level2_abundance_summer")

EFH_cols <- c(
  "2" = "#4B0082", #value <= 2
  "3" = "#006FA5", #value <= 3
  "4" = "#3CB371", #value <= 4
  "5" = "#FFFF00"  #value <= 5
)

# plot EFH map for the species
ggplot() +
  tidyterra::geom_spatvector(data = adult_alaskaplaice, ggplot2::aes(fill = as.factor(layer)), col = NA) +
  scale_fill_manual(
    name = "EFH Designation", # Optional Title for the Legend
    # map the colors to the 'layer' values (integers 1, 2, 3, 4, 5)
    values = EFH_cols,
    # define the precise text labels for the legend
    labels = c(
      "2" = "95% EFH Area (All shaded areas)",
      "3" = "75% Principal EFH Area",
      "4" = "50% Core EFH Area",
      "5" = "25% EFH Hot Spots"
    ),
    # 2. Control the legend breaks and order
    breaks = c("2", "3", "4", "5")
  )

# adjust the EFH map 
adult_alaskaplaice <- project(adult_alaskaplaice, "epsg:4326") # project EFH map onto same crs as ROMS 
sf_akplaice <- st_as_sf(adult_alaskaplaice, crs = 4326) # make sf object
sf_akplaice <- sf_akplaice |> filter(layer != "1") # remove non-EFH areas
valid_shape <- st_make_valid(sf_akplaice) # make geometry valid (not sure why)
sf_use_s2(FALSE) # don't use spherical geometry (not sure why)
valid_shape <- valid_shape |>
  st_shift_longitude() # Convert to 0-360° longitude to match ROMS data
valid_shape$layer <- as.character(valid_shape$layer) # fix EFH layer class to text instead of numbers

# plot the EFH map
ggplot() +
  geom_sf(data = valid_shape, aes(color = layer, geometry = geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_manual(
    values = EFH_cols
  ) +
  theme_bw()

# bring in temperature anomaly
temp_anoms <- data_anomaly |> st_as_sf(coords = c("lon_rho", "lat_rho")) |>
  mutate(
    anomaly_bins = case_when(anomaly >= -5 & anomaly < -2 ~"very high",
                             anomaly >= -2 & anomaly < -1.5 ~"high",
                             anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                             anomaly >= -0.5 & anomaly < 0.5~"low",
                             anomaly >= 0.5 & anomaly < 1.5~"moderate",
                             anomaly >= 1.5 & anomaly < 2~"high",
                             anomaly >=2 & anomaly <= 5~"very high")
  )
st_crs(temp_anoms)= 4326
temp_anoms = st_transform(temp_anoms, crs = st_crs(valid_shape)) # match crs of ROMS points and EFH polygons ??

# find intersects between points (ROMS) and polygons (EFH)
temp_anoms$EFH <- apply(st_intersects(valid_shape, temp_anoms, sparse = FALSE), 2, 
                                function(col) {valid_shape[which(col), ]$layer}) # not exactly sure how this code works, but it seems to have done what I want it to

# remove every point that lies outside of the species distribution
exposure <- temp_anoms |> filter(EFH == "2" | EFH == "3" | EFH == "4" | EFH == "5")

plot <- ggplot() + # plot exposure map
  geom_sf(data = exposure, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_gradientn(
    colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
    values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
    limits = c(-3, 3),
    name = "Anomaly"
  ) +
  theme_bw()
print(plot)

