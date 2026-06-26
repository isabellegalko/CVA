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
- Distributions and densities estimated from bottom trawl survey data using methods modified from Barnes et al. (2018) and Barnes et al. (2022). Standardized survey data (1990-2025) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (AFSC), NOAA and downloaded from: https://www.fisheries.noaa.gov/foss/.
- Diet-derived estimates generated using correlative spatial models and data from the Resource Ecology and Ecosystem Modeling program (Gerson et al. In prep).

**Exposure factors**
- ROMS ocean model outputs for the GOA region (temperature, salinity, phytoplankton concentration, zooplankton concentration)
- GFDL ESM4 outputs (pH, oxygen concentration, air temperature, precipitation)

## Current workflow scripts

`exposure_EFH.R` – Calculates exposure and creates associated plots for all species for which EFH maps are available and all exposure factors.

`exposure_functions.R` – Contains custom functions used in exposure_EFH.R.

`load_gfdl_data.R` – Pulls data from ESGF using OPeNDAP, separates into surface and bottom variables (as needed), and saves locally as parquet files.

`sensitivity_analysis.R` – Loads sensitivity scores from CVA workshops, calculates weighted averages, and saves in format to be used in vulnerability calculation. Also processes and computes distributional potential and directional effect scores.

## Older workflow scripts

`exposure_example.R` – Early example of calculating exposure with two species (AK plaice, walleye pollock) and two exposure factors (SST, BT).

`mapping_template.R` – Processes and creates EFH plots in the GOA using EFH data.

`temp_anomaly.R` – Processes ROMS data, computes anomalies, and creates anomaly plots for specified exposure factors.
