# Functions to create exposure maps and calculate exposure scores for the GOA CVA.
# Designed for all-species and exposure factor mapping loops.

# Default color palette for EFH designations used across maps.
EFH_cols <- c( 
  "2" = "#4B0082", #value <= 2
  "3" = "#006FA5", #value <= 3
  "4" = "#3CB371", #value <= 4
  "5" = "#FFFF00"  #value <= 5
)

sf_plot_theme = function() {
  theme_bw() +
  theme(legend.position = c(0.6, 0.17), legend.direction = "horizontal", 
        legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
        legend.key.size = unit(0.7, "cm"), legend.key.spacing = unit(0.12, "cm"),
        legend.frame = element_rect(color = "black", linewidth = 0.25),
        legend.background = element_rect(fill = "transparent"),
        legend.ticks = element_line(color = "black", linewidth = 0.25),
        strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        plot.background = element_rect(fill = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid = element_blank())
}

#' For each exposure factor sourced from ROMS, load bias-corrected parquet, 
#' filter to correct future and reference time periods, and collect data. 
#'
#' @param depth "surface" or "bottom".
#' @param variable_name Exposure factor options include: "temp", "salt", "PhL", 
#' "PhS", "Cop", "NCa", "Eup", "MZL", and "MZS"
#'
#' @assign Two data frames with ROMS data in the GOA, one for future projections 
#' (ssp585), and one for historical data (hindcasts). 
load_roms <- function(depth, variable_name) {
  variable <- open_dataset(here(paste("data/processed/all_scenarios_bias_corrected_", variable_name, "_", depth, ".parquet", sep =""))) # load projections
  future <- variable |> filter(run == "ssp585") |> # future projections
    filter(date > as.Date("2030-01-01")) |>  # Restrict dates to later than 2030
    filter(date < as.Date("2059-12-31")) |>  # Restrict to earlier than 2059
    collect()  # now load the filtered data into RAM
  hindcast <- variable |> filter(run == "hindcast") |> # hindcasts
    filter(date > as.Date("2005-01-01")) |>  # Restrict dates to later than
    filter(date < as.Date("2020-12-31")) |>  # Restrict to 2020
    collect()  # now load the filtered data into RAM
  
  assign(paste(variable_name, "_ssp585_", depth, sep = ""), future, envir = .GlobalEnv)
  assign(paste(variable_name, "_hindcast_", depth, sep = ""), hindcast, envir = .GlobalEnv)
}

#' For a particular exposure factor from ROMS (sst, bt, etc.), calculate anomaly
#' for each point in the GOA. All should be bias-corrected.
#'
#' @param future_data A data frame with variable values for future dates.
#' @param future_data A data frame with variable values for historical dates.
#'
#' @return A data frame with anomalies for each point in the GOA.
create_anomaly_roms <- function(future_data, hindcast_data) {
  # calculate climate anomaly (future mean - hindcast mean / hindcast standard deviation)
  calculate_anomaly <- future_data |> 
    filter(month == c("7", "8", "9")) |> # can filter by winter or summer months by changing month numbers
    summarize(average_future = mean(value_dc), variance_future = var(value_dc), .by = c(cell_id, lon_rho, lat_rho)) |> # use bias-corrected values: value_dc
    left_join(hindcast_data |> 
        filter(month == c("7", "8", "9")) |> 
        summarize(average_hist = mean(value), sd_hist = sd(value), variance_hist = var(value), .by = c(cell_id, lon_rho, lat_rho)), # for hindcast, value = value_dc
      by = join_by(lon_rho, lat_rho)
    ) |>
    mutate(anomaly = (average_future-average_hist)/sd_hist,
           variability = variance_hist/variance_future,
           future_change = average_future-average_hist) 
  
  anomalies <- calculate_anomaly |> 
    st_as_sf(coords = c("lon_rho", "lat_rho"))
  st_crs(anomalies)= 4326
  anomalies <- anomalies |>
  mutate(anomaly_bins = case_when(anomaly >= -5 & anomaly < -2 ~"very high",
                                   anomaly >= -2 & anomaly < -1.5 ~"high",
                                   anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                                   anomaly >= -0.5 & anomaly < 0.5 ~"low",
                                   anomaly >= 0.5 & anomaly < 1.5 ~"moderate",
                                   anomaly >= 1.5 & anomaly < 2 ~"high",
                                   anomaly >=2 & anomaly <= 5 ~"very high"))
  return(anomalies)
}

#' For a particular exposure factor from GFDL ESM (pH, o2, etc.), filter to 
#' future and reference time periods, and collect data. Construct geometry around
#' each data point in the GOA and then calculate anomaly for each cell.
#'
#' @param exposure_factor Character exposure factor name.
#' @param future_data A data frame with variable values for future dates.
#' @param future_data A data frame with variable values for historical dates.
#'
#' @return A data frame with anomalies for each point in the GOA.
create_anomaly_gfdl <- function(exposure_factor, future_data, historical_data) {
  # filter years
  future <- future_data |> filter(year >= 2030 & year <= 2059) 
  hist <- historical_data |> filter(year >= 2005 & year <= 2020) 
    
  # calculate anomaly (future mean - historical mean / historical standard deviation)
  calculate_anomaly <- future |> 
    summarize(average_future = mean(kriged), .by = c(geometry)) |>
    left_join(hist |> 
        summarize(average_hist = mean(kriged), sd_hist = sd(kriged), .by = c(geometry)),
        by = join_by(geometry)) |>
    mutate(future_change = average_future-average_hist,
           anomaly = (average_future-average_hist)/sd_hist) 
  
  spatial_data <- calculate_anomaly |>
    sf::st_as_sf(crs = 4326)
  
  anomalies <- spatial_data |>
    mutate( # set scoring categories
      anomaly_bins = case_when(anomaly >= -10 & anomaly < -2 ~"very high",
                               anomaly >= -2 & anomaly < -1.5 ~"high",
                               anomaly >= -1.5 & anomaly < -0.5 ~"moderate",
                               anomaly >= -0.5 & anomaly < 0.5 ~"low",
                               anomaly >= 0.5 & anomaly < 1.5 ~"moderate",
                               anomaly >= 1.5 & anomaly < 2 ~"high",
                               anomaly >=2 & anomaly <= 10 ~"very high"))
  return(anomalies)
}

#' For a particular exposure factor, create a plot of climate anomalies across the GOA
#' and save in plots/Anomaly plots.
#'
#' @param data A data frame with anomaly values and corresponding latitude/longitude or geometry columns.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return NA
create_anomaly_plot <- function(data, exposure_factor_name) {
  plot1 <- ggplot() + # plot climate anomaly
    geom_sf(data = data, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_gradientn(
      rescaler = function (...) {
        scales::rescale_mid(..., mid = 0)
      },
      colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
      values = scales::rescale(c(-2, -1.5, -0.5, 0.5, 1.5, 2)),
      name = "Anomaly"
    ) +
    labs(x = "Longitude",
         y = "Latitude") +
    sf_plot_theme() +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1.5))
   ggsave(paste(exposure_factor_name, "anomaly.png", sep="_"), path = here("plots/Anomaly plots"), device = "png", dpi = 300, width=8, height=5)
}

anomaly_plot_theme = function(){
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.key.size = unit(0.3, "cm"), legend.key.spacing = unit(0.1, "cm"),
        plot.margin = margin(r = 20, unit = "pt"),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        plot.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"))
}

#' For a particular exposure factor, create a series of plots:
#' 1. Average historical means across the GOA
#' 2. Historical standard deviation across the GOA
#' 3. Future change (historical - future)
#' 4. Anomalies across the GOA
#' Save plots in plots/Anomaly plots.
#'
#' @param data A data frame with anomaly values and corresponding latitude/longitude or geometry columns.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return NA
create_all_anomaly_plots <- function(data, exposure_factor_name){
  average_hist_plot <- ggplot() + 
      geom_sf(data = data, aes(color = average_hist, geometry=geometry), size = 0.5, alpha = 0.8) +
      geom_sf(data = GOA, size=0.2, fill="gray85") +
      geom_sf(data = canada, size=0.2, fill="gray95") +
      scale_x_continuous(n.breaks = 4, expand = c(0,0)) +
      scale_y_continuous(breaks = seq(52, 60, by = 4), expand = c(0,0)) +
      scale_color_gradientn(colors = c("darkslateblue",
                                       "dodgerblue1",
                                       "green3", 
                                       "yellow1",
                                       "chocolate1",  
                                       "firebrick3"),
                            name = ""
      ) +
      labs(x = "Longitude",
           y = "Latitude") +
    sf_plot_theme() +
    anomaly_plot_theme() +
      guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1.5))

  sd_hist_plot <- ggplot() + 
    geom_sf(data = data, aes(color = sd_hist, geometry=geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    scale_x_continuous(n.breaks = 4, expand = c(0,0)) +
    scale_y_continuous(breaks = seq(52, 60, by = 4), expand = c(0,0)) +
    scale_color_gradientn(colors = c("gold", "orange2", "chocolate", "firebrick4"), name = ""
    ) +
    labs(x = "Longitude") +
    sf_plot_theme() +
    anomaly_plot_theme() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1.5))
  
  future_change_plot <- ggplot() + 
    geom_sf(data = data, aes(color = future_change, geometry=geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    scale_x_continuous(n.breaks = 4, expand = c(0,0)) +
    scale_y_continuous(breaks = seq(52, 60, by = 4), expand = c(0,0)) +
    scale_color_gradientn(colors = c("darkslateblue",
                                     "dodgerblue1",
                                     "green3", 
                                     "yellow1",
                                     "chocolate1",  
                                     "firebrick3"),
                          name = ""
    ) +
    labs(x = "Longitude") +
    sf_plot_theme() +
    anomaly_plot_theme() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1.5))
  
  anomaly_plot <- ggplot() + 
    geom_sf(data = data, aes(color = anomaly, geometry=geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    scale_x_continuous(n.breaks = 4, expand = c(0,0)) +
    scale_y_continuous(breaks = seq(52, 60, by = 4), expand = c(0,0)) +
    scale_color_gradientn(
      rescaler = function (...) {
        scales::rescale_mid(..., mid = 0)
      },
      colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
      values = scales::rescale(c(-2, -1.5, -0.5, 0.5, 1.5, 2)),
      name = ""
    ) +
    labs(x = "Longitude") +
    sf_plot_theme() +
    anomaly_plot_theme() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1.5))
  
  plot_exposure <- patchwork::wrap_plots(average_hist_plot, sd_hist_plot, future_change_plot, anomaly_plot, ncol = 4)
  ggsave(filename = paste(exposure_factor_name, "_all_plots.png", sep = ""), path = here("plots/Anomaly plots"), plot = plot_exposure, device = "png", width = 8, height = 3, bg = "transparent", dpi = 400)
}

#' Prepare EFH polygons for overlap work.
#' 
#' Projects an EFH `SpatVector` to a target CRS, converts to `sf`, drops
#' non-EFH layer code "1", validates geometry, and shifts longitudes to match
#' 0-360 grids used by ROMS-style products.
#' 
#' @param path Path to EFH geodatabase files.
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param level EFH level (either 1 (presence-absence) or 2 (population percentiles)).
#'
#' @return An `sf` polygon object with valid geometry and character `layer` values.
create_EFH_layer <- function(path, species_layer, level){
  #load a specific layer
  filtered <- vect(get(path), layer = species_layer)  
  sf_use_s2(FALSE) # don't use spherical geometry 
  # adjust the EFH map 
  if(level == "2"){ # population percentiles EFH maps (level 2)
    filtered <- project(filtered, "epsg:4326") # project EFH map onto same crs as ROMS 
    sf_filtered <- st_as_sf(filtered, crs = 4326) |> # convert terra vector -> sf for sf-based spatial joins/intersections
      filter(layer == "5" | layer == "4")  |> # filter to core habitat
      st_make_valid() |> # repair invalid polygon geometries so spatial operations do not fail
      st_shift_longitude() # shift longitudes from -180..180 to 0..360 to match ROMS grid convention
    sf_filtered$layer <- as.character(sf_filtered$layer) # keep EFH designations as character labels for filtering/legend mapping
  }
  
  else{ # presence-absence EFH maps (level 1)
    filtered <- project(filtered, "epsg:4326") # project EFH map onto same crs as ROMS 
    sf_filtered <- st_as_sf(filtered, crs = 4326) |> # make sf object
      st_make_valid() |> # make geometry valid
      st_shift_longitude() }
  return(sf_filtered)
}

#' Prepare bottom trawl survey SDMs for exposure analysis.
#' 
#' Loads an .rda data structure, converts to `sf` and to a target CRS, and  
#' filters to "core" habitat area (population percentiles <50%).
#' 
#' @param path Path to .rda SDM predictions.
#'
#' @return An `sf` polygon object with valid geometry and character `layer` values.
create_BTS_SDM_layer <- function(path, species_layer){
  return(create_BTS_SDM_layer_cached(path, species_layer, TRUE))
}

# caches processed BTS data to disk and returns if already cached rather than processing again
create_BTS_SDM_layer_cached <- function(path, species_layer, use_cache = TRUE) {
  cache_path <- paste("data/filtered_sdms/", species_layer, sep = "")
  # check if this request is already cached
  if(file.exists(cache_path) && use_cache) {
    # if it is and we're using the cache, load the cached data and return it
    load(cache_path)
    return(sf_filtered)}
  
  # otherwise, load the unfiltered data and process it
  load(paste(get(path), species_layer, sep = ""))
  sf_filtered <- st_as_sf(EFH.data, coords = c("lon", "lat"), 
                          crs = 4326) |>
    st_shift_longitude() |> 
    mutate(layer = recode(EFH.area,   "25" = "5",
                          "50" = "4",
                          "75" = "3",
                          "95" = "2")) |>
    filter(layer == "5" | layer == "4") |>
    dplyr::select(!c(depth, BT, area.swept, BT.plot, fit, se.fit, pop.perc, EFH.area))
  
  # save the filtered data to the cache for next time
  save(sf_filtered, file=cache_path)
  return(sf_filtered)}

#' Prepare diet-derived SDMs for exposure analysis.
#' 
#' Loads an .rda data structure, converts to `sf` and to a target CRS, and  
#' filters to "core" habitat area (population percentiles <50%).
#' 
#' @param path Path to .rda predictions derived from diet data.
#'
#' @return An `sf` polygon object with valid geometry and character `layer` values.
create_diet_derived_layer <- function(path, species_layer){
  return(create_diet_derived_layer_cached(path, species_layer, TRUE))
}

# caches processed BTS data to disk and returns if already cached rather than processing again
create_diet_derived_layer_cached <- function(path, species_layer, use_cache = TRUE) {
  cache_path <- paste("data/filtered_sdms/", species_layer, sep = "")
  
  # check if this request is already cached
  if(file.exists(cache_path) && use_cache) {
    # if it is and we're using the cache, load the cached data and return it
    load(cache_path)
    return(sf_filtered)}
  
  # otherwise, load the unfiltered data and process it
  load(paste(get(path), species_layer, sep = ""))
  sf_filtered <- st_as_sf(prediction_data, coords = c("lon", "lat"), 
                          crs = 4326) |>
    st_shift_longitude() |> 
    mutate(layer = recode(EFH.perc,   "25" = "5",
                          "50" = "4",
                          "75" = "3",
                          "95" = "2")) |>
    filter(layer == "5" | layer == "4") |>
    dplyr::select(!c(NMFS_AREA, species, var, fit, se.fit, pop.perc, EFH.perc))
  
  # save the filtered data to the cache for next time
  save(sf_filtered, file=cache_path)
  return(sf_filtered)}

#' Prepare SDMs estimated from depth and temperature ranges for exposure analysis.
#' 
#' Loads an .rda data structure, converts to `sf` and to a target CRS.
#' 
#' @param path Path to .rda predictions derived from diet data.
#'
#' @return An `sf` polygon object with valid geometry and character `layer` values.
create_depth_temp_layer <- function(path, species_layer){
  return(create_depth_temp_layer_cached(path, species_layer, TRUE))
}

# caches processed BTS data to disk and returns if already cached rather than processing again
create_depth_temp_layer_cached <- function(path, species_layer, use_cache = TRUE) {
  cache_path <- paste("data/filtered_sdms/", species_layer, sep = "")
  
  # check if this request is already cached
  if(file.exists(cache_path) && use_cache) {
    # if it is and we're using the cache, load the cached data and return it
    load(cache_path)
    return(sf_filtered)}
  
  # otherwise, load the unfiltered data and process it
  prediction_data <- load(paste(get(path), species_layer, sep = ""))
  sf_filtered <- st_as_sf(get(prediction_data), coords = c("lon", "lat"), 
                          crs = 4326) |>
    st_shift_longitude() |> 
    dplyr::select(!c(depth, BT, species, area.swept)) |>
    mutate(layer = "5") # NOTE: this is NOT EFH layer 5, but I am making it 5 for now to work with later steps in the analysis
  
  # save the filtered data to the cache for next time
  save(sf_filtered, file=cache_path)
  return(sf_filtered)}

#' Save plot of EFH area for a single species to plots/EFH plots.
#' 
#' @param path Path to EFH or SDM files.
#' @param species_layer Either: 
#'     1. A `terra::SpatVector` containing polygons and an EFH `layer` field, or;
#'     2. An rda data structure containing points and an EFH `layer` field.
#' @param level EFH level (either 1 (presence-absence) or 2 (population percentiles)).
#'
#' @return NA
plot_species_distribution <- function(path, species_layer, level, species_name){
  # load correct species distribution
  if(path == "gdb_path" | path == "scallop_path"){
    sf_filtered = create_EFH_layer(path, species_layer, level)}
  else if(path == "bts_sdm_path"){
    sf_filtered = create_BTS_SDM_layer(path, species_layer)}
  else if(path == "diet_derived_path"){
    sf_filtered = create_diet_derived_layer(path, species_layer)}
  else if(path == "depth_temp_path"){
    sf_filtered = create_depth_temp_layer(path, species_layer)}
  
  #create plot
  if(path == "gdb_path"){ # for EFH predictions
  plot <- ggplot() +
    geom_sf(data = sf_filtered, aes(color = layer, fill = layer, geometry = geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_manual(values = EFH_cols,
                       guide = "none") +
    scale_fill_manual(values = EFH_cols,
                      labels = c(
                        "2" = "95% EFH Area",
                        "3" = "75% Principal EFH Area",
                        "4" = "50% Core EFH Area",
                        "5" = "25% EFH Hot Spots"),
                      # 2. Control the legend breaks and order
                      breaks = c("2", "3", "4", "5")) +
    labs(x = "Longitude",
         y = "Latitude",
         fill = "EFH area") +
    sf_plot_theme() +
    theme(panel.grid = element_blank())
  }
  else if(path == "scallop_path" | path == "depth_temp_path"){ # for presence-absence SDMs
  plot <- ggplot() +
    geom_sf(data = sf_filtered, aes(geometry = geometry), size = 0.5, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    labs(x = "Longitude",
         y = "Latitude") +
    sf_plot_theme()
  }
  else if(path == "bts_sdm_path" | path == "diet_derived_path"){ # for EFH predictions from BTS or diet data
  plot <- ggplot() +
      geom_sf(data = sf_filtered, aes(geometry = geometry, color = layer), size = 0.3, alpha = 0.8) +
      geom_sf(data = GOA, size=0.2, fill="gray85") +
      geom_sf(data = canada, size=0.2, fill="gray95") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      scale_color_manual(
        values = c("4" = "#3CB371",
                   "5" = "#FFFF00"),
        labels = c(
          "4" = "50% Core EFH Area",
          "5" = "25% EFH Hot Spots"),
        name = "Population \n Percentile",
        breaks = c("4", "5")) +
      labs(x = "Longitude",
           y = "Latitude") +
      sf_plot_theme()
  }
  ggsave(paste(species_name, "distribution.png", sep="_"), device = "png", path = here("plots/EFH plots"), plot = plot, width=8, height=5, dpi = 300)
}

#' Keep points that intersect species EFH polygons.
#' 
#' Reprojects anomaly points to EFH CRS, computes point-in-polygon overlaps,
#' and stores matching a EFH layer value per point.
#' 
#' @param path Path to EFH or SDM files.
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param level EFH level (either 1 (presence-absence) or 2 (population percentiles)).
#' @param species_name Character species name.
#' @param anomaly_data A data frame with calculated anomalies for each point in the GOA.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return Filtered `sf` points representing exposure within EFH.
create_overlap <- function(path, species_layer, level, species_name, anomaly_data, exposure_factor_name){
  if(path == "gdb_path" | path == "scallop_path"){
    sf_filtered = create_EFH_layer(path, species_layer, level)

    anomaly_data = st_transform(anomaly_data, crs = st_crs(sf_filtered)) # transform points to exact CRS used by EFH polygons
    
    # spatial join - assign EFH layer to each ROMS point
    # st_intersects() creates a polygon x point logical matrix:
    # apply(..., 2, ...) loops across columns (one point at a time)
    # for each point-column, which(col) finds intersecting polygon row indices,
    # then sf_filtered[which(col), ]$layer returns that point's EFH layer label(s)
    
    # for species with EFH designations
    if(level == "2"){ 
      anomaly_data$EFH <- apply(st_intersects(sf_filtered, anomaly_data, sparse = FALSE), 2, 
                                function(col) {sf_filtered[which(col), ]$layer}) 
      exposure <- anomaly_data |> 
        filter(EFH == "4" | EFH == "5") # filter to core habitat (50% EFH Area)
    }
    # for species with EFH presence-absence maps (e.g., weathervane scallop)
    else{ 
      anomaly_data$EFH <- apply(st_intersects(sf_filtered, anomaly_data, sparse = FALSE), 2, 
                                function(col) {sf_filtered[which(col), ]$geometry}) 
      
      exposure <- anomaly_data |> mutate(na=map_lgl(.x = EFH, .f = is_empty)) # identify rows not not included in the presence-absence polygon
      exposure <- exposure |> filter(na == FALSE)    
    }
  }
  else if(path == "bts_sdm_path"){ # for SDMs
    sf_filtered = create_BTS_SDM_layer(path, species_layer)
    nearest_points <- st_nearest_feature(sf_filtered, anomaly_data)
    exposure <- cbind(sf_filtered, st_drop_geometry(anomaly_data)[nearest_points, ])
  }
  else if(path == "diet_derived_path"){ # for SDMs
    sf_filtered = create_diet_derived_layer(path, species_layer)
    nearest_points <- st_nearest_feature(sf_filtered, anomaly_data)
    exposure <- cbind(sf_filtered, st_drop_geometry(anomaly_data)[nearest_points, ])
  }
  else if(path == "depth_temp_path"){ # for SDMs
    sf_filtered = create_depth_temp_layer(path, species_layer)
    nearest_points <- st_nearest_feature(sf_filtered, anomaly_data)
    exposure <- cbind(sf_filtered, st_drop_geometry(anomaly_data)[nearest_points, ])
  }
  return(exposure)
}

#' Bin anomalies into categories from low to very high.
#' 
#' Assigns each point an exposure level: low, moderate, high, or very high 
#' according to previous CVA methods (see Loughran et al. 2025).
#' 
#' @param path Path to EFH or SDM files.
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param level EFH level (either 1 (presence-absence) or 2 (population percentiles)).
#' @param species_name Character species name.
#' @param anomaly_data A data frame with calculated anomalies for each point in the GOA.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return Data frame with count of anomalies in each exposure category.
assign_exposure_levels <- function(data){
  exposure <- data
  # assign scoring categories from low - very high to the anomaly values
  exposure_plot <- exposure |> 
    mutate(
      exposure_score = ifelse(anomaly >= -0.5 & anomaly <= 0.5, "low", ifelse((anomaly < -0.5 & anomaly >= -1.5) | (anomaly > 0.5 & anomaly <= 1.5), "moderate", ifelse((anomaly < -1.5 & anomaly >= -2) | (anomaly > 1.5 & anomaly <= 2), "high", "very_high")))) |>
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
      prop = count / total)
  return(exposure_plot)
}

#' Calculate an exposure score for a particular species and exposure factor.
#' 
#' @param path Path to EFH or SDM files.
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param level EFH level (either 1 (presence-absence) or 2 (population percentiles)).
#' @param species_name Character species name.
#' @param anomaly_data A data frame with calculated anomalies for each point in the GOA.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return An integer (an exposure score for a single species/exposure factor).
calculate_exposure_score <- function(original_exposure_data, path, species_layer, level, species_name, anomaly_data, exposure_factor_name){
  exposure_plot <- assign_exposure_levels(original_exposure_data)
  # calculate weighted mean
  exp_fact_mean <- exposure_plot |>
    dplyr::select(!c(prop,total)) |>
    pivot_wider(names_from = "exposure_score",
                values_from = "count") |>
    mutate(weighted_mean = ((low*1) + (moderate*2) + (high*3) + (very_high*4))/(low + moderate + high + very_high))
  exp_fact_mean$weighted_mean <- round(exp_fact_mean$weighted_mean, digits = 2)
  
  return(exp_fact_mean$weighted_mean)
}

#' Bin anomalies into categories from low to very high (according to CVA methods)
#' and save three plots in plots/Exposure plots/species name:
#' 1. An exposure map of anomalies across the GOA within a particular species' EFH area.
#' 2. A histogram of binned anomalies from the exposure map.
#' 3. The distribution of the exposure scores across the four scoring categories.
#' 
#' @param path Path to EFH or SDM files.
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param level EFH level (either 1 (presence-absence) or 2 (population percentiles)).
#' @param species_name Character species name.
#' @param anomaly_data A data frame with calculated anomalies for each point in the GOA.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return NA
create_exposure_plots <- function(path, species_layer, level, species_name, anomaly_data, exposure_factor_name) {
  # only do first time!
  # dir.create(here("plots/Exposure plots/"), recursive = TRUE)
  
  original_exposure_data <- create_overlap(path, species_layer, level, species_name, anomaly_data, exposure_factor_name)
  exposure_plot <- assign_exposure_levels(original_exposure_data)
  
  # plot exposure map
  plot2 <- ggplot() + 
    geom_sf(data = original_exposure_data, aes(color = anomaly, geometry=geometry), size = 1, alpha = 0.8) +
    geom_sf(data = GOA, size=0.2, fill="gray85") +
    geom_sf(data = canada, size=0.2, fill="gray95") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_gradientn(
      rescaler = function (...) {
        scales::rescale_mid(..., mid = 0)
      },
      colors = c("purple", "blue", "cyan", "green", "yellow", "orange", "red"), # set colors for scoring categories
      values = scales::rescale(c(-2, -1.5, -0.5, 0.5, 1.5, 2)),
      name = "") +
    labs(x = "Longitude",
         y = "Latitude",
         color = "Anomaly") +
    theme_bw() +
    theme(panel.grid = element_blank(),
      legend.position = c(0.6, 0.15), legend.direction = "horizontal", 
      legend.text = element_text(size = 8), legend.title = element_text(size = 10), 
      legend.key.size = unit(0.7, "cm"), legend.key.spacing = unit(0.12, "cm"),
      legend.frame = element_rect(color = "black", linewidth = 0.25),
      legend.background = element_rect(fill = "transparent"),
      legend.ticks = element_line(color = "black", linewidth = 0.25),
      strip.text = element_text(hjust = 0, size = 10),
      strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
      plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  ggsave(paste(species_name, exposure_factor_name, "exposure_map.png", sep="_"), path = here("plots/Exposure plots/"), plot = plot2, device = "png", width = 8, height = 5, dpi = 300)
  
  # make histogram of anomalies
  plot3 <- ggplot(original_exposure_data) +
    geom_histogram(aes(x = anomaly, y = after_stat(count / sum(count)), fill = anomaly_bins), binwidth = 0.25, boundary = 0, linewidth = 0.25, colour="black", show.legend = FALSE) +
    scale_fill_manual(values = c("low" = "green", 
                                 "moderate" = "yellow", 
                                 "high" = "orange", 
                                 "very high" = "red")) +
    scale_x_continuous(
      breaks = seq(-10, 10, by = 1), 
      limits = c(-3, 3),
      expand = c(0,0)
    ) +
    xlab("Anomaly") +
    ylab("Proportion") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  ggsave(paste(species_name, exposure_factor_name, "anomaly_histogram.png", sep="_"), path = here("plots/Exposure plots/"), plot = plot3, device = "png")
  
  exp_fact_mean <- calculate_exposure_score(original_exposure_data, path, species_layer, level, species_name, anomaly_data, exposure_factor_name)
  
  # plot distribution of exposure scores
  plot4 <- ggplot(exposure_plot) +
    geom_col(mapping = aes(x = exposure_score, y = prop, fill = exposure_score), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
    scale_x_discrete(labels = c("L", "M", "H", "V")) +
    ylab("Proportion") +
    xlab("Exposure Score") +
    scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
    annotate("text", x = I(0.8), y = I(0.9), label = paste("Exposure = ", exp_fact_mean, sep = "")) + # add exposure score in top right corner
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  ggsave(paste(species_name, exposure_factor_name, "exposure_scores.png", sep="_"), path = here("plots/Exposure plots/"), plot = plot4, device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)
  
}

# plot distribution of exposure scores for all exposure factors for each species
exposure_histogram_series <- function(path, species_layer, EFH_level, species_name){
  exposure_plot_list = list() # create list
  # select correct exposure factors for each species
  exposure_factors_list <- layer_names[i,6:13]
  these_exposure_factors <- as.list(as.data.frame(t(exposure_factors_list)))
  these_exposure_factors <- lapply(these_exposure_factors, function(x) x[!is.na(x)])
  
  # create data frame to put exposure score numbers in
  species_scores <- data.frame(
    exposure_factor = unlist(these_exposure_factors$V1)
  )
  
  for(i in 1:length(these_exposure_factors$V1)){
    original_exposure_data <- create_overlap(path, species_layer, EFH_level, species_name, anomaly[[paste(these_exposure_factors$V1[i], "_anomaly", sep = "")]], these_exposure_factors$V1[i])
    exposure_plot <- assign_exposure_levels(original_exposure_data)
    exposure_plot <- exposure_plot |>
      mutate(exposure_factor = these_exposure_factors$V1[i])
    # assign(paste(exposure_factors[i], "_exposure_plot", sep = ""), exposure_plot, envir = .GlobalEnv)
    exposure_plot_list[[i]] <- exposure_plot # add it to the list
    
    # calculate exposure score for each exposure factor
    score <- calculate_exposure_score(original_exposure_data, path, species_layer, EFH_level, species_name, anomaly[[paste(these_exposure_factors$V1[i], "_anomaly", sep = "")]], these_exposure_factors$V1[i])
    species_scores[i,"score"] <- score
  }
  # create single data frame with every df from the for loop
  # combine all exposure factor results into a single df per species
  group_exposure_plot <- do.call(rbind, exposure_plot_list) 
  
  group_plot <- ggplot(group_exposure_plot) +
    geom_col(mapping = aes(x = exposure_score, y = prop, fill = exposure_score), 
             position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
    geom_text(data = species_scores, aes(x = I(0.9), y = I(0.92), label = score), size = 3) +
    facet_wrap(~exposure_factor) +
    scale_x_discrete(labels = c("L", "M", "H", "V")) +
    ylab("Proportion") +
    xlab("Exposure Score") +
    scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
    #annotate("text", x = I(0.8), y = I(0.8), label = paste("Exposure = ", exp_fact_mean, sep = "")) + # add exposure score in top right corner
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          plot.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  
  ggsave(paste(species_name, "all_exposure_scores.png", sep="_"), path = here("plots/Exposure plots/"), plot = group_plot, device = "png", width = 7, height = 5, bg = "transparent", dpi = 300)
}
