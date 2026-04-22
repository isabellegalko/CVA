# ============================================================================
# EFH Mapping Script - Version 2 (All Species)
# ============================================================================
# Purpose: Create EFH + overlap maps for all species and write:
#   1) Per-species plot files in plots/species-name/
#   2) One combined multi-page PDF (3 maps per species on one page)
# Author: Mason Smith, Isabelle Galko, Chelsey Beese

# Load required packages -------------------------------------------------------
library(here)
library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(tibble)

# Source helper functions ------------------------------------------------------
source(here::here("map_helper_functions.r"))

# Project path helper (works with or without {here}) --------------------------
project_path <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    return(here::here(...))
  }
  file.path(getwd(), ...)
}

# Inputs ----------------------------------------------------------------------
gdb_path <- here::here("data", "GOA_groundfish_2023.gdb")
plots_base_dir <- here::here("plots")
combined_pdf <- here::here("plots", "all-species-maps.pdf")

# Pattern used to identify adult EFH layers in the geodatabase
layer_pattern <- "^GOA_adult_.*_efh_level2_abundance_summer$"

# Discover EFH layers and species names ---------------------------------------
layers <- terra::vector_layers(gdb_path)
species_layers <- layers[grepl(layer_pattern, layers)]

if (length(species_layers) == 0L) {
  stop("No EFH layers matched pattern: ", layer_pattern)
}

species_names <- species_layers |>
  gsub("^GOA_adult_", "", x = _) |>
  gsub("_efh_level2_abundance_summer$", "", x = _) |>
  gsub("_", " ", x = _)

# Build anomaly points once (shared by all species) ---------------------------
temp_anoms_sf <- prepare_temp_anoms_sf(data_anomaly)

# Main species loop ------------------------------------------------------------
species_pages <- vector("list", length(species_layers))
names(species_pages) <- species_names

species_summary <- vector("list", length(species_layers))

for (i in seq_along(species_layers)) {
  layer_name <- species_layers[[i]]
  species_name <- species_names[[i]]

  # Read and prepare species EFH polygons
  efh_vect <- terra::vect(gdb_path, layer = layer_name)
  efh_sf <- prepare_efh_sf(efh_vect)

  # Spatial overlap: anomaly points inside EFH
  temp_with_efh <- add_efh_overlap(temp_anoms_sf, efh_sf, efh_col_name = "EFH")
  exposure <- filter_exposure_points(temp_with_efh, efh_col_name = "EFH")

  # Build the three species plots
  p_efh <- plot_EFH(efh_sf)
  p_check <- check_EFH(efh_sf, coast)
  p_exposure <- plot_exposure(exposure, coast)

  # Save individual species plots to plots/species-name/
  species_files <- save_species_plots(
    species_name = species_name,
    plot_list = list(efh = p_efh, check = p_check, exposure = p_exposure),
    base_dir = plots_base_dir
  )

  # Build one-page (3 panel) species summary for combined PDF
  species_pages[[i]] <- compose_species_page(
    species_name = species_name,
    efh_plot = p_efh,
    check_plot = p_check,
    exposure_plot = p_exposure
  )

  species_summary[[i]] <- tibble::tibble(
    species = species_name,
    layer = layer_name,
    n_exposure_points = nrow(exposure),
    plot_dir = species_plot_dir(species_name, base_dir = plots_base_dir),
    efh_plot_file = species_files[[1]],
    check_plot_file = species_files[[2]],
    exposure_plot_file = species_files[[3]]
  )
}

# Write one giant all-species PDF (one page per species) ----------------------
write_species_summary_pdf(
  page_plots = species_pages,
  output_file = combined_pdf,
  width = 16,
  height = 6
)

# Final object returned by script ---------------------------------------------
run_summary <- dplyr::bind_rows(species_summary)
run_summary