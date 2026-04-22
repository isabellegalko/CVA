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

# ...existing code...

# adjust the EFH map so it can be compared with ROMS point data
adult_alaskaplaice <- project(adult_alaskaplaice, "epsg:4326") # reproject EFH polygons to WGS84 (same base CRS as ROMS)
sf_akplaice <- st_as_sf(adult_alaskaplaice, crs = 4326) |> # convert terra vector -> sf for sf-based spatial joins/intersections
  filter(layer != "1") |> # drop layer 1 (non-EFH area)
  st_make_valid() # repair invalid polygon geometries so spatial operations do not fail
sf_use_s2(FALSE) # use planar geometry engine instead of spherical s2 for these operations
sf_akplaice <- sf_akplaice |> st_shift_longitude() # shift longitudes from -180..180 to 0..360 to match ROMS grid convention
sf_akplaice$layer <- as.character(sf_akplaice$layer) # keep EFH designations as character labels for filtering/legend mapping

# quick check plot: EFH polygons over coastline
# color corresponds to EFH layer (2, 3, 4, 5)
ggplot() +
  geom_sf(data = sf_akplaice, aes(color = layer, geometry = geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_manual(
    values = EFH_cols
  ) +
  theme_bw()

# bring in temperature anomaly points from ROMS output
# convert lon/lat columns to an sf point object, then classify anomaly magnitude into bins
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
st_crs(temp_anoms)= 4326 # explicitly assign CRS to ROMS points (WGS84)
temp_anoms = st_transform(temp_anoms, crs = st_crs(sf_akplaice)) # transform points to exact CRS used by EFH polygons

# spatial join logic: assign EFH layer to each ROMS point
# st_intersects(sf_akplaice, temp_anoms, sparse = FALSE) creates a polygon x point logical matrix:
# - rows = EFH polygons
# - columns = ROMS points
# - TRUE means the point is inside/intersects that polygon
# apply(..., 2, ...) loops across columns (one point at a time)
# for each point-column, which(col) finds intersecting polygon row indices,
# then sf_akplaice[which(col), ]$layer returns that point's EFH layer label(s)
temp_anoms$EFH <- apply(st_intersects(sf_akplaice, temp_anoms, sparse = FALSE), 2, 
                                function(col) {sf_akplaice[which(col), ]$layer})

# keep only ROMS points that fall within EFH layers 2-5 (inside species distribution)
exposure <- temp_anoms |> filter(EFH == "2" | EFH == "3" | EFH == "4" | EFH == "5")

# final exposure map:
# - points colored by continuous anomaly value
# - coastline added for geographic reference
# - custom diverging-like color scale centered around near-zero anomalies
plot <- ggplot() +
  geom_sf(data = exposure, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
  geom_sf(data = coast, color = "black", linewidth = 0.3) +
  scale_color_gradientn(
    colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # anomaly color ramp from cold to warm extremes
    values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)), # where each color sits along the -3 to 3 range
    limits = c(-3, 3),
    name = "Anomaly"
  ) +
  theme_bw()
print(plot)


