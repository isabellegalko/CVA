# Gulf of Alaska Climate Vulnerability Assessment

This repository contains code and analysis for the Gulf of Alaska (GOA) Climate Vulnerability Assessment (CVA).

## Definitions
- **Climate vulnerability** – how a species will be impacted by climate change, composed of sensitivity and exposure
- **Exposure** – conditions that a species experiences in its environment, quantified through the overlap of species distributions and environmental projections
- **Exposure factors** – environmental variables that are of interest to the list of species
- **Sensitivity** – the intrinsic susceptability of a species to change
- **Sensitivity attributes** – a set of life history characteristics that characterize a species' potential response to climate change

## Data sources 
**Species distributions** 
- Essential Fish Habitat (EFH) predictions for federally-managed species in the GOA. File geodatabases for species in the Groundfish and Scallop Fishery Management Plans (FMP) were downloaded from: https://www.fisheries.noaa.gov/resource/map/alaska-essential-fish-habitat-efh-species-shapefiles.
- `data/bts_sdms`: Distributions and densities estimated from bottom trawl survey data using methods modified from Barnes et al. (2018) and Barnes et al. (2022). Standardized survey data (1990-2025) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (AFSC), NOAA and downloaded from: https://www.fisheries.noaa.gov/foss/.
- `data/diet_sdms`: Diet-derived estimates generated using correlative spatial models and data from the Resource Ecology and Ecosystem Modeling program (Gerson et al. In prep).
- `data/depth_temp_sdms`: Data-poor distribution estimates using depth and temperature ranges.

**Exposure factors**
- ROMS ocean model outputs: 
  - temperature, salinity, phytoplankton concentration, zooplankton concentration
- GFDL ESM4 outputs:
  - `data/pH`: pH (ocean acidification)
  - `data/o2`: oxygen concentration
  - `data/tas`: air temperature
  - `data/pr`: precipitation
  

## Current workflow scripts

1. `load_gfdl_data.R` – Pulls data from ESGF using OPeNDAP, separates into surface and bottom variables (as needed), and saves locally as parquet files.

2. `exposure_EFH.R` – Calculates exposure and creates associated plots for all species for which EFH maps are available and all exposure factors.
  - `exposure_functions.R` – Contains custom functions.

3. `Vulnerability.R` – Calculates vulnerability by combining exposure and sensitivity scores. Creates plots.
