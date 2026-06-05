This repository contains code and analysis for the Gulf of Alaska (GOA) Climate Vulnerability Assessment (CVA).

### Definitions
**Climate vulnerability** – how a species will be impacted by climate change, composed of sensitivity and exposure

**Exposure** – conditions that a species experiences in its environment, quantified through the overlap of species distributions and environmental projections

**Exposure factors** – environmental variables that are of interest to the list of species

### Data sources 
**Species distributions** – geodatabase of Essential Fish Habitat (EFH) predictions for the GOA Region

**Exposure factors**
- ROMS ocean model outputs for the GOA region (temperature, salinity, phytoplankton concentration, zooplankton concentration)
- GFDL ESM4 outputs (pH, oxygen concentration, air temperature, precipitation)

### Current workflow scripts

**exposure_EFH.R** – Calculates exposure and creates associated plots for all species for which EFH maps are available and all exposure factors.

**exposure_functions.R** – Contains custom functions used in exposure_EFH.R.

**load_gfdl_data.R** – Pulls data from ESGF using OPeNDAP and saves locally as parquet files.

### Older workflow scripts

**exposure_example.R** – Early example of calculating exposure with two species (AK plaice, walleye pollock) and two exposure factors (SST, BT).

**mapping_template.R** – Processes and creates EFH plots in the GOA using EFH data.

**temp_anomaly.R** – Processes ROMS data, computes anomalies, and creates anomaly plots for specified exposure factors.
