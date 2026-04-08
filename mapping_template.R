
#title: EFH mapping script
#author: Mason Smith
#date: 8 April 2026
#editor: visual

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

#plot
ggplot() +
  tidyterra::geom_spatvector(data = adult_alaskaplaice, ggplot2::aes(fill = as.factor(layer)), col = NA) +
  scale_fill_manual(
    name = "EFH Designation", # Optional Title for the Legend
    # map the colors to the 'layer' values (integers 1, 2, 3, 4, 5)
    values = c(
      "2" = "#4B0082", #value <= 2
      "3" = "#006FA5", #value <= 3
      "4" = "#3CB371", #value <= 4
      "5" = "#FFFF00"  #value <= 5
    ),
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

# crds(adult_alaskaplaice)

adult_alaskaplaice <- project(adult_alaskaplaice, "epsg:4326") 
sf_akplaice <- st_as_sf(adult_alaskaplaice) 
sf_akplaice <- sf_akplaice |> filter(layer != "1")

shifted_akp <- st_shift_longitude(sf_akplaice)
valid_shape <- st_make_valid(shifted_akp)
sf_use_s2(TRUE)

temp_anoms <- data_anomaly |> st_as_sf(coords = c("lon_rho", "lat_rho"), crs = 4326)

points_in_polygons <- st_filter(temp_anoms, valid_shape, .predicate = st_intersects)

# points_in_polygons <- intersect(temp_anoms, shifted_akp)
