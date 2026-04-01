#' Read ROMS Grid Data with Spatial Subset
#' 
#' @description
#' Reads the ROMS grid file and subsets it to match the indices used in the 
#' pre-processed NetCDF files (xi_rho: 67-226, eta_rho: 237-450)
#' 
#' @param grid_file Path to the ROMS grid NetCDF file (default: 'data/NEP_grid_5a.nc')
#' 
#' @return A tibble with columns: xi_rho_new, eta_rho_new, lon_rho, lat_rho, h (depth)
#' 
#' @examples
#' roms_grid <- read_roms_grid()
#' 
#' @import tidync
#' @import dplyr
#' @importFrom purrr pluck map_df
read_roms_grid <- function(grid_file = here::here('data', 'NEP_grid_5a.nc')) {
  
  # Open grid file
  roms <- tidync(grid_file)
  
  # Find available grids
  roms_vars <- hyper_grids(roms) %>%
    pluck("grid") %>%
    purrr::map_df(function(x) {
      roms %>% activate(x) %>% hyper_vars() %>% 
        mutate(grd = x)
    })
  
  # Get rho grid
  latlon_rhogrd <- roms_vars %>% filter(name == "lat_rho") %>% pluck('grd')
  roms_rho <- roms %>% 
    activate(latlon_rhogrd) %>% 
    hyper_tibble() %>%
    dplyr::select(lon_rho, lat_rho, xi_rho, eta_rho, h)# %>% 
    #mutate(lon_rho = lon_rho - 360) # Convert to -180 to 180
  
  # Subset to match the pre-processed data indices
  # Original indices: xi_rho 67-226, eta_rho 237-450
  roms_rho <- roms_rho %>%
    filter(between(xi_rho, 67, 226), between(eta_rho, 237, 450)) %>%
    mutate(xi_rho_new = xi_rho - 66,  # Renumber starting from 1
           eta_rho_new = eta_rho - 236) %>%
    dplyr::select(-xi_rho, -eta_rho)
  
  return(roms_rho)
}


#' Process Single Annual NetCDF File
#' 
#' @description
#' Reads one annual NetCDF file containing surface and bottom layer data for
#' multiple variables, with optional filtering by year range and depth.
#' Files outside the year range are skipped before reading.
#' 
#' @param ncfile Name of the NetCDF file to process (e.g., "annual_1990.nc")
#' @param data_dir Directory containing the NetCDF files
#' @param roms_grid ROMS grid data from read_roms_grid()
#' @param variables Character vector of variable names to extract. 
#'   Default is all 10 variables: temp, salt, PhS, PhL, MZS, MZL, Cop, NCa, Eup, Det
#' @param min_year Minimum year to include in output (required)
#' @param max_year Maximum year to include in output (required)
#' @param maxdepth Maximum depth (h) to include in meters (default: 1000)
#' 
#' @return A tibble with columns: date, variable, layer, lon_rho, lat_rho, value
#'   Returns NULL if file year is outside the specified range
#' 
#' @examples
#' grid <- read_roms_grid()
#' data <- process_annual_file("annual_1990.nc", "data/annual_files", grid, # though usually ran over batches of files
#'                            min_year = 1990, max_year = 1990)
#' data <- process_annual_file("annual_2010.nc", "data/annual_files", grid, 
#'                            min_year = 2010, max_year = 2010, maxdepth = 500)
#' 
#' @import tidync
#' @import dplyr
#' @import tidyr
#' @import ncdf4
#' @import lubridate
process_annual_file <- function(ncfile, data_dir, roms_grid, 
                                variables = c("temp", "salt", "PhS", "PhL", 
                                              "MZS", "MZL", "Cop", "NCa", 
                                              "Eup", "Det"),
                                min_year = NA,
                                max_year = NA,
                                maxdepth = 1000,
                                mask = goa_mask) {
  
  
  # Check if min_year and max_year are provided
  if (is.na(min_year) || is.na(max_year)) {
    stop("Both min_year and max_year must be provided.")
  }
  
  # Extract year from filename (e.g., "annual_1990.nc" -> 1990)
  file_year <- as.numeric(gsub("annual_(\\d{4})\\.nc", "\\1", ncfile))
  
  # Check if file year is within range - skip if not
  if (file_year < min_year || file_year > max_year) {
    print(paste("Skipping file:", ncfile, "(outside year range)"))
    return(NULL)
  }
  
  print(paste("Processing file:", ncfile))
  
  # Full path to file
  filepath <- here::here(data_dir, ncfile)
  
  # Open with tidync
  nc <- tidync(filepath)
  
  # Get all variables in long format
  nc_data <- nc %>% 
    hyper_tibble(na.rm = FALSE) %>%
    dplyr::select(xi_rho, eta_rho, ocean_time, s_rho, all_of(variables))
  
  # Convert to long format for easier plotting
  nc_data_long <- nc_data %>%
    pivot_longer(cols = all_of(variables), 
                 names_to = "variable", 
                 values_to = "value") %>%
    drop_na(value)  # Remove NA values (land cells)
  
  # rename xi and eta columns to reflect the fact that these are no longer the original indices but they have been reset to 1
  nc_data_long <- nc_data_long %>% rename(xi_rho_new = xi_rho, eta_rho_new = eta_rho)
  
  # Join with grid to get coordinates
  nc_data_long <- nc_data_long %>%
    left_join(roms_grid, by = c("xi_rho_new", "eta_rho_new"))
  
  # Filter by depth using maxdepth parameter
  nc_data_long <- nc_data_long %>% filter(h < maxdepth)
  
  # eliminate xi and eta points outside the ROMS mask
  # this slows down the function but it will produce much smaller masks
  nc_data_long <- nc_data_long %>%
    mutate(idx_drop = paste(xi_rho_new, eta_rho_new, sep = "_")) %>%
    filter(!idx_drop %in% mask) %>%
    dplyr::select(-idx_drop)
  
  # Convert ocean_time to dates
  nc_file <- nc_open(filepath)
  time_data <- ncvar_get(nc_file, "ocean_time")
  time_units <- ncatt_get(nc_file, "ocean_time", "units")$value
  time_parts <- strsplit(time_units, " ")[[1]]
  ref_date_str <- paste(time_parts[3:length(time_parts)], collapse = " ")
  nc_close(nc_file)
  
  # Create date lookup
  dates <- data.frame(
    ocean_time = time_data,
    date = as.POSIXct(time_data, origin = ref_date_str, tz = "UTC")
  ) %>%
    mutate(date = as.Date(date))
  
  # Add dates to data
  nc_data_long <- nc_data_long %>%
    left_join(dates, by = "ocean_time")
  
  # Add layer labels
  nc_data_long <- nc_data_long %>%
    mutate(layer = case_when(
      s_rho == max(s_rho) ~ "surface",
      s_rho == min(s_rho) ~ "bottom",
      TRUE ~ "other"
    ))
  
  # drop unneeded cols and transform to factor where possible
  nc_data_long <- nc_data_long %>%
    dplyr::select(date, variable, layer, lon_rho, lat_rho, value) %>%
    mutate(
      variable = factor(variable),
      layer = factor(layer)
    )
  
  return(nc_data_long)
}

#' Create Parquet File from ROMS Annual Data
#'
#' This function processes ROMS annual data files with monthly time steps
#' for specified runs and years, combining them into a single parquet file.
#'
#' @param run Character string specifying the run type. Must be one of:
#'   \itemize{
#'     \item "hindcast" - Historical hindcast run (default years: 1990-2020)
#'     \item "historical" - Historical run (default years: 1980-2014)
#'     \item "ssp585" - SSP5-8.5 scenario (default years: 2015-2099)
#'     \item "ssp245" - SSP2-4.5 scenario (default years: 2015-2099)
#'     \item "ssp126" - SSP1-2.6 scenario (default years: 2015-2099)
#'   }
#' @param min_year Integer. Minimum year to process. If NULL (default), uses
#'   run-specific defaults: hindcast (1990), historical (1980), 
#'   ssp585/ssp245/ssp126 (2015).
#' @param max_year Integer. Maximum year to process. If NULL (default), uses
#'   run-specific defaults: hindcast (2020), historical (2014),
#'   ssp585/ssp245/ssp126 (2099).
#' @param variables Character vector of variables to process. Default includes
#'   all available variables: "temp" (temperature), "salt" (salinity), 
#'   "PhS" (small phytoplankton), "PhL" (large phytoplankton), 
#'   "MZS" (small microzooplankton), "MZL" (large microzooplankton),
#'   "Cop" (copepods), "NCa" (neocalanus), "Eup" (euphausiids), 
#'   "Det" (detritus).
#' @param maxdepth Numeric. Maximum depth to process. Default is 1000.
#' @param mask Character vector. Mask for which ROMS points to drop. Default is
#'   "goa_mask", which will be read from \code{idx_to_drop.txt}.
#'
#' @return Invisibly returns the path to the saved parquet file.
#' 
#' @details
#' The function reads NetCDF files from the directory structure 
#' \code{data/annual_files/{run}/}, processes them using \code{process_annual_file()},
#' and saves the combined output to \code{data/processed/}. If an output file 
#' already exists, a new file will be created with a numbered suffix (e.g., 
#' "hindcast_annual_data (1).parquet").
#' 
#' The function requires:
#' \itemize{
#'   \item A \code{functions.R} file with the \code{process_annual_file()} and 
#'         \code{read_roms_grid()} functions
#'   \item An \code{idx_to_drop.txt} file containing the ROMS grid mask
#'   \item NetCDF files in the appropriate \code{data/annual_files/{run}/} directory
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process hindcast with all defaults
#' create_parquet_file(run = "hindcast")
#' 
#' # Process historical with custom years
#' create_parquet_file(run = "historical", min_year = 1985, max_year = 2010)
#' 
#' # Process SSP5-8.5 scenario with subset of variables
#' create_parquet_file(run = "ssp585", variables = c("temp", "salt"))
#' 
#' # Process with custom depth limit
#' create_parquet_file(run = "hindcast", maxdepth = 500)
#' 
#' # Process with all custom parameters
#' create_parquet_file(
#'   run = "ssp585",
#'   min_year = 2050,
#'   max_year = 2080,
#'   variables = c("temp", "PhS", "PhL", "Cop"),
#'   maxdepth = 750
#' )
#' 
#' }
create_parquet_file <- function(run,
                                min_year = NULL,
                                max_year = NULL,
                                variables = c("temp", "salt", "PhS", "PhL", 
                                              "MZS", "MZL", "Cop", "NCa", 
                                              "Eup", "Det"),
                                maxdepth = 1000,
                                mask = "goa_mask") {
  
  # Validate run argument
  valid_runs <- c("hindcast", "historical", "ssp585", "ssp245", "ssp126")
  if (!run %in% valid_runs) {
    stop("Invalid 'run' argument. Must be one of: ", 
         paste(valid_runs, collapse = ", "))
  }
  
  # Set default years based on run if not provided
  if (is.null(min_year)) {
    min_year <- switch(run,
                       hindcast = 1990,
                       historical = 1980,
                       ssp585 = 2015,
                       ssp245 = 2015,
                       ssp126 = 2015)
  }
  
  if (is.null(max_year)) {
    max_year <- switch(run,
                       hindcast = 2020,
                       historical = 2014,
                       ssp585 = 2099,
                       ssp245 = 2099,
                       ssp126 = 2099)
  }
  
  # Generate base output filename
  base_outfile <- paste0("data/processed/", run, "_annual_data.parquet")
  
  # Handle file name collision
  outfile <- base_outfile
  counter <- 1
  while (file.exists(outfile)) {
    outfile <- paste0("data/processed/", run, "_annual_data (", counter, ").parquet")
    counter <- counter + 1
  }
  
  if (outfile != base_outfile) {
    message(paste("Output file already exists. Writing to:", basename(outfile)))
  }
  
  # Source functions
  source("functions.R")
  
  # Read ROMS grid
  grid <- read_roms_grid()
  
  # Read in mask for which ROMS point to drop
  if (mask == "goa_mask") {
    goa_mask <- scan("idx_to_drop.txt", "character", sep = " ")
  } else {
    goa_mask <- mask
  }
  
  # Read all annual files and bind them together
  data_dir <- paste0("data/annual_files/", run)
  nc_files <- list.files(data_dir, 
                         pattern = "annual_.*\\.nc$", 
                         full.names = FALSE)
  
  if (length(nc_files) == 0) {
    stop("No NetCDF files found in ", data_dir)
  }
  
  message(paste("Found", length(nc_files), "files to process"))
  message(paste("Processing run:", run))
  message(paste("Year range:", min_year, "to", max_year))
  message(paste("Variables:", paste(variables, collapse = ", ")))
  
  # Initialize empty list to store results
  all_data_list <- list()
  
  # Loop through files
  for (i in seq_along(nc_files)) {
    file <- nc_files[i]
    
    # Process file
    data <- process_annual_file(
      ncfile = file,
      data_dir = data_dir,
      roms_grid = grid,
      min_year = min_year,
      max_year = max_year,
      variables = variables,
      maxdepth = maxdepth,
      mask = goa_mask
    )
    
    # Add to list
    all_data_list[[i]] <- data
    
    # Progress message
    message(paste("Completed", i, "of", length(nc_files)))
    
    # Optional: garbage collection every few files to free up memory
    if (i %% 5 == 0) {
      gc()
    }
  }
  
  # Bind all data together
  all_data <- bind_rows(all_data_list)
  rm(all_data_list)
  gc()
  
  # Sort by date
  all_data <- all_data %>%
    arrange(date, variable, layer) %>%
    mutate(run = factor(run))
  
  # Summary
  message(paste("Total rows:", nrow(all_data)))
  message(paste("Date range:", min(all_data$date), "to", max(all_data$date)))
  message(paste("Variables:", paste(unique(all_data$variable), collapse = ", ")))
  
  # Save the combined dataset
  write_parquet(all_data, outfile)
  message(paste("Data saved to:", outfile))
  
  # Return the output file path invisibly
  invisible(outfile)
}


#' Create Spatial Maps for Surface and Bottom Layers
#' 
#' @description
#' Creates faceted spatial maps showing surface and bottom layers for a 
#' specific variable, year, and month with coastline overlay
#' 
#' @param data Data frame from reading annual files
#' @param variable Variable name to plot (e.g., "temp", "salt")
#' @param year Year to plot (e.g., 2010)
#' @param month Month to plot (1-12)
#' @param coastline sf object containing coastline data from rnaturalearth
#' @param title Plot title (optional, auto-generated if NULL)
#' 
#' @return A ggplot object with faceted maps for surface and bottom layers
#' 
#' @examples
#' # Load coastline data once
#' library(rnaturalearth)
#' coast <- ne_coastline(scale = "medium", returnclass = "sf") %>%
#'   st_crop(xmin = -170, xmax = -130, ymin = 50, ymax = 62) %>%
#'   st_shift_longitude()
#' 
#' # Create map for June 2010
#' plot_spatial_map(all_data, "temp", 2010, 6, coast)
#' 
#' # Create map for different variable
#' plot_spatial_map(all_data, "salt", 2015, 8, coast)
#' 
#' @import ggplot2
#' @import dplyr
#' @import sf
#' @import rnaturalearth
plot_spatial_map <- function(data, variable, year, month, coastline, title = NULL, psize = 0.5) {
  
  title <- paste0(variable, " in ", year, "-", month)
  
  dat <- data %>%
    mutate(yr = year(date),
           mo = month(date)) %>%
    filter(variable == !!variable, yr == year, mo == month)
  
  # Create plot
  p <- dat %>%
    st_as_sf(coords = c("lon_rho", "lat_rho"), crs = 4326) %>%
    ggplot() +
    geom_sf(aes(color = value), size = psize, alpha = 0.8) +
    geom_sf(data = coastline, color = "black", linewidth = 0.3) +
    facet_grid(run~layer) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = title,
         x = "Longitude",
         y = "Latitude",
         color = variable) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Create Time Series Plot with Ribbon
#' 
#' @description
#' Creates a faceted time series plot showing spatial mean with ribbons 
#' representing spatial variability (5th-95th percentiles) for surface 
#' and bottom layers
#' 
#' @param data Data frame from reading annual files
#' @param variable Variable name to plot (e.g., "temp", "salt")
#' @param start_year Starting year for the plot (optional, e.g., 2000)
#' @param end_year Ending year for the plot (optional, e.g., 2010)
#' @param title Plot title (optional, auto-generated if NULL)
#' 
#' @return A ggplot object with faceted panels for surface and bottom layers
#' 
#' @examples
#' # Full time series
#' plot_time_series(all_data, "temp")
#' 
#' # Specific year range
#' plot_time_series(all_data, "temp", start_year = 2000, end_year = 2010)
#' 
#' # Just years after 2005
#' plot_time_series(all_data, "salt", start_year = 2005)
#' 
#' @import ggplot2
#' @import dplyr
#' @import lubridate
plot_time_series <- function(data, variable, start_year = NULL, end_year = NULL, title = NULL) {
  
  # Filter by variable
  plot_data <- data %>%
    filter(variable == !!variable)
  
  # Filter by year range if provided
  if (!is.null(start_year)) {
    start_date <- as.Date(paste0(start_year, "-01-01"))
    plot_data <- plot_data %>%
      filter(date >= start_date)
  }
  
  if (!is.null(end_year)) {
    end_date <- as.Date(paste0(end_year, "-12-31"))
    plot_data <- plot_data %>%
      filter(date <= end_date)
  }
  
  # Calculate spatial statistics with 5th and 95th percentiles
  plot_data <- plot_data %>%
    group_by(date, layer, run) %>%
    summarise(mean_value = mean(value, na.rm = TRUE),
              lower = quantile(value, 0.05, na.rm = TRUE),
              upper = quantile(value, 0.95, na.rm = TRUE),
              .groups = "drop")
  
  # Set title
  if (is.null(title)) {
    year_range <- if (!is.null(start_year) || !is.null(end_year)) {
      paste0(" (", 
             ifelse(is.null(start_year), "", start_year),
             "-",
             ifelse(is.null(end_year), "", end_year),
             ")")
    } else {
      ""
    }
    title <- paste0("Spatial mean ", variable, " over time", year_range, 
                    "\n(shaded area: 5th-95th percentiles)")
  }
  
  # Create plot with facets
  p <- ggplot(plot_data, aes(x = date, y = mean_value, color = run, fill = run)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line(linewidth = 0.8) +
    facet_grid(~layer) +
    labs(title = title,
         x = "Date",
         y = paste("Mean", variable)) +
    theme_bw()
  
  return(p)
}

#' Bias Correct ROMS Projection Data Using Delta Method
#'
#' @description
#' Applies delta correction to ROMS projection data based on the difference
#' between hindcast (data-assimilated) and historical (free-running) during
#' a reference period. Correction is applied at the cell-month level to
#' preserve spatial and seasonal patterns in the bias.
#'
#' @param hindcast Data frame with hindcast data (1990-2020)
#' @param historical Data frame with historical data (1980-2014)
#' @param projection Data frame with projection data (2015-2099)
#' @param ref_years Vector of reference years for calculating bias. Default is 1991:2014
#' @param lognormal Logical. Apply correction on log scale? Default FALSE
#' @param use_sd Logical. Use ratio of standard deviations in correction? Default FALSE
#' @param include_hindcast Logical. Splice hindcast into output time series? Default TRUE
#'
#' @return Data frame with bias-corrected time series containing:
#'   - Bias-corrected historical (1980 through Jan 1990)
#'   - Hindcast (Feb 1990 through Dec 2020) if include_hindcast = TRUE
#'   - Bias-corrected projection (Jan 2021 through 2099)
#'   Columns: date, variable, layer, lon_rho, lat_rho, value, value_dc, run, source_run
#'
#' @details
#' The function calculates monthly mean temperatures for each grid cell during
#' the reference period (1991-2014) for both hindcast and historical runs.
#' These are used as correction factors.
#'
#' Delta correction formula (without SD scaling):
#' value_corrected = mean_hindcast + (value - mean_historical)
#'
#' With SD scaling:
#' value_corrected = mean_hindcast + (sd_hindcast/sd_historical) * (value - mean_historical)
#'
#' For lognormal variables, the correction is applied on the log scale.
#'
#' The output includes:
#' - value: original value
#' - value_dc: delta-corrected value
#' - run: the projection scenario name (e.g., "ssp585")
#' - source_run: which run the data originally came from ("historical", "hindcast", "projection")
#'
#' @examples
#' \dontrun{
#' # Load data
#' hindcast <- open_dataset("data/processed/hindcast_annual_data.parquet") %>%
#'   filter(layer == "surface", variable == "temp") %>%
#'   collect()
#'
#' historical <- open_dataset("data/processed/historical_annual_data.parquet") %>%
#'   filter(layer == "surface", variable == "temp") %>%
#'   collect()
#'
#' ssp585 <- open_dataset("data/processed/ssp585_annual_data.parquet") %>%
#'   filter(layer == "surface", variable == "temp") %>%
#'   collect()
#'
#' # Apply bias correction
#' corrected_ssp585 <- bias_correct_roms(
#'   hindcast = hindcast,
#'   historical = historical,
#'   projection = ssp585,
#'   use_sd = FALSE,
#'   include_hindcast = TRUE
#' )
#' }
#'
#' @export
bias_correct_roms <- function(
    hindcast,
    historical,
    projection,
    ref_years = 1991:2014,
    lognormal = FALSE,
    use_sd = FALSE,
    include_hindcast = TRUE) {
  
  # Add temporal variables if not present
  if(!"year" %in% names(hindcast)) {
    hindcast <- hindcast %>% mutate(year = year(date), month = month(date))
  }
  if(!"year" %in% names(historical)) {
    historical <- historical %>% mutate(year = year(date), month = month(date))
  }
  if(!"year" %in% names(projection)) {
    projection <- projection %>% mutate(year = year(date), month = month(date))
  }
  
  # Create cell identifiers for joining
  hindcast <- hindcast %>%
    mutate(cell_id = paste(round(lon_rho, 3), round(lat_rho, 3), sep = "_"))
  
  historical <- historical %>%
    mutate(cell_id = paste(round(lon_rho, 3), round(lat_rho, 3), sep = "_"))
  
  projection <- projection %>%
    mutate(cell_id = paste(round(lon_rho, 3), round(lat_rho, 3), sep = "_"))
  
  # Store original projection run name
  proj_run_name <- unique(projection$run)[1]
  
  # Convert to log if needed
  if(lognormal) {
    projection <- projection %>% mutate(value = log(value))
    historical <- historical %>% mutate(value = log(value))
    hindcast <- hindcast %>% mutate(value = log(value))
  }
  
  # Calculate mean and SD for each cell-month-variable-layer during reference period
  # - Historical
  correction_hist <- historical %>%
    filter(year %in% ref_years) %>%
    group_by(cell_id, lon_rho, lat_rho, month, variable, layer) %>%
    summarise(
      mean_hist = mean(value, na.rm = TRUE),
      sd_hist = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # - Hindcast
  correction_hind <- hindcast %>%
    filter(year %in% ref_years) %>%
    group_by(cell_id, lon_rho, lat_rho, month, variable, layer) %>%
    summarise(
      mean_hind = mean(value, na.rm = TRUE),
      sd_hind = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge correction factors
  correction_factors <- correction_hist %>%
    left_join(correction_hind, by = c("cell_id", "lon_rho", "lat_rho", "month", "variable", "layer"))
  
  # Apply correction to projection
  projection_corrected <- projection %>%
    left_join(correction_factors, by = c("cell_id", "lon_rho", "lat_rho", "month", "variable", "layer"))
  
  # Apply correction to historical
  historical_corrected <- historical %>%
    left_join(correction_factors, by = c("cell_id", "lon_rho", "lat_rho", "month", "variable", "layer"))
  
  # Calculate corrected values
  if(!lognormal) {
    if(use_sd) {
      # With SD scaling
      projection_corrected <- projection_corrected %>%
        mutate(value_dc = mean_hind + (sd_hind / sd_hist * (value - mean_hist)))
      
      historical_corrected <- historical_corrected %>%
        mutate(value_dc = mean_hind + (sd_hind / sd_hist * (value - mean_hist)))
    } else {
      # Without SD scaling (default)
      projection_corrected <- projection_corrected %>%
        mutate(value_dc = mean_hind + (value - mean_hist))
      
      historical_corrected <- historical_corrected %>%
        mutate(value_dc = mean_hind + (value - mean_hist))
    }
  } else {
    # Lognormal correction
    if(use_sd) {
      projection_corrected <- projection_corrected %>%
        mutate(value_dc = exp(mean_hind + (sd_hind / sd_hist * (value - mean_hist))))
      
      historical_corrected <- historical_corrected %>%
        mutate(value_dc = exp(mean_hind + (sd_hind / sd_hist * (value - mean_hist))))
    } else {
      projection_corrected <- projection_corrected %>%
        mutate(value_dc = exp(mean_hind + (value - mean_hist)))
      
      historical_corrected <- historical_corrected %>%
        mutate(value_dc = exp(mean_hind + (value - mean_hist)))
    }
    
    # Back-transform original values
    projection_corrected <- projection_corrected %>% mutate(value = exp(value))
    historical_corrected <- historical_corrected %>% mutate(value = exp(value))
    hindcast <- hindcast %>% mutate(value = exp(value))
  }
  
  # Clean up correction factor columns
  projection_corrected <- projection_corrected %>%
    select(date, year, month, variable, layer, lon_rho, lat_rho, cell_id, value, value_dc, run) %>%
    mutate(source_run = proj_run_name)
  
  historical_corrected <- historical_corrected %>%
    select(date, year, month, variable, layer, lon_rho, lat_rho, cell_id, value, value_dc, run) %>%
    mutate(source_run = "historical")
  
  # Prepare hindcast (value_dc = value for hindcast)
  hindcast_prepared <- hindcast %>%
    mutate(
      value_dc = value,
      source_run = "hindcast"
    ) %>%
    select(date, year, month, variable, layer, lon_rho, lat_rho, cell_id, value, value_dc, run, source_run)
  
  # Assemble final time series
  if(include_hindcast) {
    # Historical (<Feb 1990) + Hindcast (Feb 1990 - Dec 2020) + Projection (>2020)
    full_time_series <- bind_rows(
      historical_corrected %>% filter(date < as.Date("1990-02-01")),
      hindcast_prepared,
      projection_corrected %>% filter(year > max(hindcast_prepared$year))
    ) %>%
      arrange(date) %>%
      mutate(run = proj_run_name)  # Label entire series with projection scenario name
    
  } else {
    # Historical + Projection (no hindcast splice)
    full_time_series <- bind_rows(
      historical_corrected,
      projection_corrected
    ) %>%
      arrange(date) %>%
      mutate(run = proj_run_name)
  }
  
  return(full_time_series)
}
