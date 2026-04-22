# Functions to create EFH and exposure maps for CVA workflows.
# Designed for all-species mapping loops.

# Default color palette for EFH designations used across maps.
default_efh_cols <- c(
  "2" = "#4B0082",
  "3" = "#006FA5",
  "4" = "#3CB371",
  "5" = "#FFFF00"
)

#' Prepare EFH polygons for overlap work
#'
#' Projects an EFH `SpatVector` to a target CRS, converts to `sf`, drops
#' non-EFH layer code "1", validates geometry, and shifts longitudes to match
#' 0-360 grids used by ROMS-style products.
#'
#' @param efh_vect A `terra::SpatVector` containing EFH polygons and a `layer` field.
#' @param target_crs EPSG code used to project EFH polygons (default `4326`).
#' @param drop_layer Character layer value to remove (default `"1"`).
#' @param shift_longitude Logical; if `TRUE`, applies `sf::st_shift_longitude()`.
#'
#' @return An `sf` polygon object with valid geometry and character `layer` values.
prepare_efh_sf <- function(
    efh_vect,
    target_crs = 4326,
    drop_layer = "1",
    shift_longitude = TRUE
) {
  sf::sf_use_s2(FALSE)

  efh_sf <- terra::project(efh_vect, paste0("epsg:", target_crs)) |>
    sf::st_as_sf(crs = target_crs) |>
    dplyr::mutate(layer = as.character(layer)) |>
    dplyr::filter(layer != drop_layer) |>
    sf::st_make_valid()

  if (shift_longitude) {
    efh_sf <- efh_sf |> sf::st_shift_longitude()
  }

  efh_sf
}

#' Convert ROMS anomaly table to sf points
#'
#' Builds point geometry from longitude/latitude columns and adds a categorical
#' anomaly bin for quick severity mapping.
#'
#' @param data_anomaly A data frame with anomaly values and coordinate columns.
#' @param lon_col Name of longitude column in `data_anomaly` (default `"lon_rho"`).
#' @param lat_col Name of latitude column in `data_anomaly` (default `"lat_rho"`).
#' @param crs CRS EPSG code assigned to output points (default `4326`).
#'
#' @return An `sf` point object with an `anomaly_bins` column.
prepare_temp_anoms_sf <- function(data_anomaly, lon_col = "lon_rho", lat_col = "lat_rho", crs = 4326) {
  data_anomaly |>
    sf::st_as_sf(coords = c(lon_col, lat_col), crs = crs) |>
    dplyr::mutate(
      anomaly_bins = dplyr::case_when(
        anomaly >= -5 & anomaly < -2 ~ "very high",
        anomaly >= -2 & anomaly < -1.5 ~ "high",
        anomaly >= -1.5 & anomaly < -0.5 ~ "moderate",
        anomaly >= -0.5 & anomaly < 0.5 ~ "low",
        anomaly >= 0.5 & anomaly < 1.5 ~ "moderate",
        anomaly >= 1.5 & anomaly < 2 ~ "high",
        anomaly >= 2 & anomaly <= 5 ~ "very high",
        TRUE ~ NA_character_
      )
    )
}

#' Add EFH overlap classification to anomaly points
#'
#' Reprojects anomaly points to EFH CRS, computes point-in-polygon overlaps,
#' and stores matching EFH layer value(s) per point.
#'
#' @param temp_anoms_sf `sf` points of ROMS anomalies.
#' @param efh_sf `sf` polygons with a `layer` column.
#' @param efh_col_name Output column name for overlap labels (default `"EFH"`).
#'
#' @return Input `sf` points with an added EFH overlap column.
add_efh_overlap <- function(temp_anoms_sf, efh_sf, efh_col_name = "EFH") {
  temp_anoms_sf <- sf::st_transform(temp_anoms_sf, sf::st_crs(efh_sf))

  idx <- sf::st_intersects(temp_anoms_sf, efh_sf, sparse = TRUE)

  temp_anoms_sf[[efh_col_name]] <- vapply(
    idx,
    function(i) {
      if (length(i) == 0L) {
        NA_character_
      } else {
        paste(efh_sf$layer[i], collapse = ";")
      }
    },
    FUN.VALUE = character(1)
  )

  temp_anoms_sf
}

#' Keep points that intersect species EFH polygons
#'
#' Filters anomaly points to retain rows where EFH overlap exists.
#'
#' @param temp_anoms_with_efh `sf` points with an EFH overlap column.
#' @param efh_col_name Column name storing EFH overlap labels.
#'
#' @return Filtered `sf` points representing exposure within EFH.
filter_exposure_points <- function(temp_anoms_with_efh, efh_col_name = "EFH") {
  temp_anoms_with_efh |>
    dplyr::filter(!is.na(rlang::.data[[efh_col_name]]))
}

#' Plot EFH designation polygons
#'
#' Generates an EFH polygon map with manual legend labels/colors for layers 2-5.
#'
#' @param efh_data `sf` or `terra::SpatVector` polygon data with `layer` field.
#' @param efh_cols Named vector of colors keyed by EFH layer code.
#'
#' @return A `ggplot` object.
plot_EFH <- function(efh_data, efh_cols = default_efh_cols) {
  base_plot <- ggplot2::ggplot()

  if (inherits(efh_data, "sf")) {
    base_plot <- base_plot +
      ggplot2::geom_sf(data = efh_data, ggplot2::aes(fill = as.factor(layer)), color = NA)
  } else {
    base_plot <- base_plot +
      tidyterra::geom_spatvector(data = efh_data, ggplot2::aes(fill = as.factor(layer)), col = NA)
  }

  base_plot +
    ggplot2::scale_fill_manual(
      name = "EFH Designation",
      values = efh_cols,
      labels = c(
        "2" = "95% EFH Area (All shaded areas)",
        "3" = "75% Principal EFH Area",
        "4" = "50% Core EFH Area",
        "5" = "25% EFH Hot Spots"
      ),
      breaks = c("2", "3", "4", "5")
    ) +
    ggplot2::theme_bw()
}

#' Quick EFH geometry check against coastlines
#'
#' Creates a diagnostic map of EFH polygon outlines over coastline geometry.
#'
#' @param sf_data Species EFH polygons as `sf`.
#' @param coast_data Coastline geometry as `sf`.
#' @param efh_cols Named vector of colors keyed by EFH layer code.
#'
#' @return A `ggplot` object.
check_EFH <- function(sf_data, coast_data, efh_cols = default_efh_cols) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf_data, ggplot2::aes(color = layer, geometry = geometry), size = 0.5, alpha = 0.8) +
    ggplot2::geom_sf(data = coast_data, color = "black", linewidth = 0.3) +
    ggplot2::scale_color_manual(values = efh_cols) +
    ggplot2::theme_bw()
}

#' Plot species exposure (anomaly values inside EFH)
#'
#' Draws anomaly points as colored `sf` points on coastline background.
#'
#' @param exposure_data `sf` points filtered to EFH overlap.
#' @param coast_data Coastline geometry as `sf`.
#'
#' @return A `ggplot` object.
plot_exposure <- function(exposure_data, coast_data) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = exposure_data, ggplot2::aes(color = anomaly, geometry = geometry), size = 0.5, alpha = 0.8) +
    ggplot2::geom_sf(data = coast_data, color = "black", linewidth = 0.3) +
    ggplot2::scale_color_gradientn(
      colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"),
      values = scales::rescale(c(-3, -2, -1.5, -0.5, 0, 0.5, 1.5, 2, 3)),
      limits = c(-3, 3),
      name = "Anomaly"
    ) +
    ggplot2::theme_bw()
}

#' Build a clean species folder path under a base plots directory
#'
#' @param species_name Character species name.
#' @param base_dir Base directory for all plot output (default `"plots"`).
#'
#' @return Character path to the species output directory.
species_plot_dir <- function(species_name, base_dir = "plots") {
  species_slug <- tolower(species_name)
  species_slug <- gsub("[^a-z0-9]+", "-", species_slug)
  species_slug <- gsub("(^-|-$)", "", species_slug)

  fs::path(base_dir, species_slug)
}

#' Save the three per-species map files into plots/species-name/
#'
#' Saves named plots as PNG files in each species folder.
#' Expected names are typically: `efh`, `check`, and `exposure`.
#'
#' @param species_name Character species name.
#' @param plot_list Named list of ggplot objects.
#' @param base_dir Base output directory for plots.
#' @param width Plot width passed to `ggsave`.
#' @param height Plot height passed to `ggsave`.
#' @param units Units passed to `ggsave`.
#' @param dpi DPI passed to `ggsave`.
#'
#' @return Character vector of saved file paths.
save_species_plots <- function(
    species_name,
    plot_list,
    base_dir = "plots",
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
) {
  out_dir <- species_plot_dir(species_name, base_dir = base_dir)
  fs::dir_create(out_dir)

  plot_names <- names(plot_list)
  if (is.null(plot_names) || any(plot_names == "")) {
    plot_names <- paste0("plot_", seq_along(plot_list))
  }

  out_files <- fs::path(out_dir, paste0(plot_names, ".png"))

  purrr::walk2(
    .x = plot_list,
    .y = out_files,
    .f = ~ ggplot2::ggsave(filename = .y, plot = .x, width = width, height = height, units = units, dpi = dpi)
  )

  out_files
}

#' Compose one single-page species summary with three maps
#'
#' Arranges three map panels (EFH, check, exposure) on one page for the
#' species-level page used in the combined all-species PDF.
#'
#' @param species_name Character species name used in page title.
#' @param efh_plot EFH map `ggplot` object.
#' @param check_plot EFH check map `ggplot` object.
#' @param exposure_plot Exposure map `ggplot` object.
#'
#' @return A patchwork plot object representing one species page.
compose_species_page <- function(species_name, efh_plot, check_plot, exposure_plot) {
  patchwork::wrap_plots(efh_plot, check_plot, exposure_plot, nrow = 1) +
    patchwork::plot_annotation(title = species_name)
}

#' Write a multi-page PDF from species page plots
#'
#' Writes a single PDF where each plot in `page_plots` is a separate page.
#'
#' @param page_plots A list of one-page ggplot/patchwork objects.
#' @param output_file Path to the combined PDF file.
#' @param width Page width (inches).
#' @param height Page height (inches).
#'
#' @return Path to the generated PDF file.
write_species_summary_pdf <- function(page_plots, output_file, width = 16, height = 6) {
  fs::dir_create(fs::path_dir(output_file))

  grDevices::pdf(file = output_file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  purrr::walk(page_plots, ~ print(.x))

  output_file
}
