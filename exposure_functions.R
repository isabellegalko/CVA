# Functions to create exposure maps and calculate exposure scores for the GOA CVA.

# colors for plotting EFH areas
EFH_cols <- c( 
  "2" = "#4B0082", #value <= 2
  "3" = "#006FA5", #value <= 3
  "4" = "#3CB371", #value <= 4
  "5" = "#FFFF00"  #value <= 5
)

#' For a particular exposure factor, create a plot of climate anomalies from ROMS data.
#'
#' @param data A data frame with anomaly values and coordinate columns.
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return A ggplot object.
create_anomaly_plot <- function(data, exposure_factor_name) {
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

#' Prepare EFH polygons for overlap work
#' 
#' Projects an EFH `SpatVector` to a target CRS, converts to `sf`, drops
#' non-EFH layer code "1", validates geometry, and shifts longitudes to match
#' 0-360 grids used by ROMS-style products.
#' 
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#'
#' @return An `sf` polygon object with valid geometry and character `layer` values.
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

#' Create plot of EFH area for a single species
#' 
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#'
#' @return A ggplot object.
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



#' Add EFH overlap classification to anomaly points
#' 
#' Reprojects anomaly points to EFH CRS, computes point-in-polygon overlaps,
#' and stores matching a EFH layer value per point.
#' 
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param species_name Character species name.
#' @param anomaly_data 
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return Data frame
create_overlap <- function(species_layer, species_name, anomaly_data, exposure_factor_name){
  sf_filtered <- create_EFH_layer(species_layer)
  
  anomaly_data = st_transform(anomaly_data, crs = st_crs(sf_filtered)) # match crs of ROMS points and EFH polygons ??
  
  # find intersects between points (ROMS) and polygons (EFH)
  anomaly_data$EFH <- apply(st_intersects(sf_filtered, anomaly_data, sparse = FALSE), 2, 
                            function(col) {sf_filtered[which(col), ]$layer}) # not exactly sure how this code works, but it seems to have done what I want it to
  
  # remove every point that lies outside of the species distribution
  exposure <- anomaly_data |> filter(EFH == "4" | EFH == "5") # filter to core habitat (50% EFH Area)
  
  return(exposure)
}

#' Bins anomalies into categories from low to very high (according to CVA methods)
#' 
#' Assigns each point an exposure level: low, moderate, high, or very high.
#' 
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param species_name Character species name.
#' @param anomaly_data 
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return 
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

#' Calculate an exposure score for a particular species and exposure factor.
#' 
#' Calculates a weighted mean of the anomalies.
#' 
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param species_name Character species name.
#' @param anomaly_data 
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return An integer.
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


#' Bins anomalies into categories from low to very high (according to CVA methods)
#' 
#' Assigns each point an exposure level: low, moderate, high, or very high.
#' 
#' @param species_layer A `terra::SpatVector` containing EFH polygons and a `layer` field. 
#' @param species_name Character species name.
#' @param anomaly_data 
#' @param exposure_factor_name Character exposure factor name.
#'
#' @return A ggplot object: an exposure map of anomalies across the GOA within a particular species' EFH area.
#' @return A ggplot object: histogram of binned anomalies from the exposure map.
#' @return A ggplot object: distribution of the exposure scores across the four scoring categories.
create_exposure_plots <- function(species_layer, species_name, anomaly_data, exposure_factor_name) {
  #dir.create(here(paste("plots/Exposure plots/", species_name, "/", exposure_factor_name, sep ="")), recursive = TRUE)
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
