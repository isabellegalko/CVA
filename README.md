# Gulf of Alaska Climate Vulnerability Assessment

This repository contains code and analysis for the Gulf of Alaska (GOA) Climate Vulnerability Assessment (CVA).

## Definitions
- **Climate vulnerability** – how a species will be impacted by climate change, composed of sensitivity and exposure
- **Exposure** – conditions that a species experiences in its environment, quantified through the overlap of species distributions and environmental projections
- **Exposure factors** – environmental variables that are of interest to the list of species
- **Sensitivity** – the intrinsic susceptibility of a species to change
- **Sensitivity attributes** – a set of life history characteristics that characterize a species' potential response to climate change

## Data  
**Species list** 
- `exposure_assignments.xlsx`: Full species list, distribution file paths, and exposure assignments.

**Species distributions** 
- Essential Fish Habitat (EFH) predictions for federally-managed species in the GOA. File geodatabases for species in the Groundfish and Scallop Fishery Management Plans (FMP) were downloaded from: https://www.fisheries.noaa.gov/resource/map/alaska-essential-fish-habitat-efh-species-shapefiles.
  - These files are too large for GitHub. Please download them yourself and include the geodatabases in the `data` folder.
- `data/bts_sdms`: Distributions and densities estimated from bottom trawl survey data using methods modified from Barnes et al. (2018) and Barnes et al. (2022). Standardized survey data (1990-2025) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (AFSC), NOAA and downloaded from: https://www.fisheries.noaa.gov/foss/.
- `data/diet_sdms`: Diet-derived estimates generated using correlative spatial models and data from the Resource Ecology and Ecosystem Modeling program (Gerson et al. In prep).
- `data/depth_temp_sdms`: Data-poor distribution estimates using depth and temperature ranges.

**Exposure factors**
- ROMS ocean model outputs: 
  - Variables: temperature, salinity, phytoplankton concentration, zooplankton concentration
  - These files are too large for GitHub. Please pull the data from ROMS yourself and include them in the `data` folder.
- GFDL ESM4 outputs:
  - `data/pH`: pH (ocean acidification)
  - `data/o2`: oxygen concentration
  - `data/tas`: air temperature
  - `data/pr`: precipitation
- Marine heatwave index: 
  - `data/mhw`: Cumulative intensity of marine heatwaves (see https://github.com/somros/MHWI_GOACLIM; Rovellini et al. in prep)
  
## Workflow

`00_functions.R` - Contains custom functions.

`01_sensitivity_analysis.R` – Processes sensitivity scoring from CVA workshops. Calculates distribution change potential, directional effect, and data quality scores for each species and functional group. Saves sensitivity results to `results`.

`02_load_gfdl_outputs.R` – Pulls files from ESGF using OPeNDAP, separates into surface and bottom variables (as needed), and saves locally as parquet files.

`03_exposure_analysis.R` – Calculates exposure scores and creates associated plots. Saves exposure results to `results`.

`04_vulnerability.R` – Calculates vulnerability by combining exposure and sensitivity scores. Synthesizes trends in sensitivity, exposure, and vulnerability across functional groups. Conducts leave-one-out analysis and bootstrapping of vulnerability scores. Creates plots.
