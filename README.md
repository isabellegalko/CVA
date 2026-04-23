This repository contains code and analysis for the Gulf of Alaska (GOA) Climate Vulnerability Assessment (CVA).

### Definitions
**Climate vulnerability** – how a species will be impacted by climate change, composed of sensitivity and exposure

**Exposure** – conditions that a species experiences in its environment, quantified through the overlap of species distributions and environmental projections

**Exposure factors** – environmental variables that are of interest to the list of species

### Data sources 
**Species distributions** – geodatabase of Essential Fish Habitat (EFH) predictions for the GOA Region

**Exposure factors** – ROMS ocean model output data for the GOA region

### Notable workflow scripts

**exposure_EFH.R** – Calculates exposure and creates associated plots for all species for which EFH maps are available and several exposure factors (SST, BT, and salinity).

**exposure_functions.R** – Contains custom functions used in exposure_EFH.R.

**exposure_example.R** – Early example of calculating exposure with two species (AK plaice, walleye pollock) and two exposure factors (SST, BT).

**mapping_template.R** – Processes and creates EFH plots in the GOA using EFH data.

**temp_anomaly.R** – Processes ROMS data, computes anomalies, and creates anomaly plots for specified exposure factors.
