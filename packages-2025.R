# ------------------------------------------------------------------------
# File: packages.R
# Purpose: Install R packages after a fresh R setup
# Author: Chelsey Beese
# Date created: 7th May 2022
# Date modified: 24th Nov 2025
# ------------------------------------------------------------------------

# Package groupings for improved workflow
# 1. Environment & Editor
env_editor <- c(
    "distill", "here", "installr",
    "languageserver", "remotes", "rstudioapi",
    "usethis"
)

# 2. Development Tools
dev_tools <- c(
    "covr", "devtools", "lintr",
    "pkgdown", "roxygen2", "testthat",
    "tryCatchLog"
)

# 3. Reporting & Publishing
reporting <- c(
    "blogdown", "bookdown", "kableExtra",
    "knitcitations", "knitr", "magick",
    "rmarkdown", "tinytex"
)

# 4. Data Wrangling & Analysis
data <- c(
    "car", "data.table", "dplyr",
    "FSA", "gmodels", "janitor",
    "lme4", "lubridate", "MASS",
    "MESS", "multcomp", "nlme",
    "readxl", "reshape2", "rfishbase",
    "tidymodels", "tidyverse", "vegan",
    "xtable", "plyr", "abind",
    "reshape"
)

# 5. Visualization
viz <- c(
    "cowplot", "fmsb", "gganimate",
    "ggbeeswarm", "ggforce", "ggmap",
    "ggpattern", "ggpubr", "ggtext",
    "ggthemes", "patchwork", "plotly",
    "plotrix", "png", "ragg", "rgl",
    "RColorBrewer", "wesanderson", "scales"
)

# 6. Special/Experimental
special <- c(
    "htmlwidgets", "shiny", "shinyBS",
    "shinyjs", "XLConnect", "xlsx"
)

# 7. Alaska work
alaska <- c("sf", "terra", "raster", "easyNCDF", "DirichletReg")


# Install CRAN packages
install.packages(env_editor)
install.packages(dev_tools)
install.packages(reporting)
install.packages(data)
install.packages(viz)
install.packages(special)
install.packages(alaska)

# 7. Project-specific or GitHub packages (sorted by name)
install.packages("mizer")
github_packages <- list(
    # ggvegan = "gavinsimpson/ggvegan",
    mizerExperimental = "sizespectrum/mizerExperimental",
    # mizerShelf = "gustavdelius/mizerShelf" # Uncomment if needed
    mizerMR = "sizespectrum/mizerMR",
    mizerReef = "cmbeese/mizerReef",
    akgfmaps = "afsc-gap-products/akgfmaps", # maps for Alaska
    coldpool = "afsc-gap-products/coldpool", # calculate the cold pool index, mean sea surface temperature, and mean bottom temperature
    akfishcondition = "afsc-gap-products/akfishcondition", # fish condition metrics for Alaska fish species
    catchfunction = "amandafaig/catchfunction" # various catch functions for size spectrum models
)

# Install GitHub packages
library(remotes)
for (pkg in github_packages) {
    remotes::install_github(pkg, build_vignettes = TRUE)
}

# ------------------------------------------------------------------------
# Package Descriptions (all packages installed above)
# ------------------------------------------------------------------------

# 1. Environment & Editor
# distill         | Create and publish scientific and technical articles
# here            | Simplifies file paths relative to the project directory
# installr        | Tools for installing and updating R and Rtools (Windows only)
# languageserver  | Language Server Protocol support for R in IDEs
# remotes         | Install R packages from remote repositories (GitHub, GitLab, etc.)
# rstudioapi      | Access RStudio IDE features programmatically
# usethis         | Automate package and project setup tasks

# 2. Development Tools
# covr            | Test coverage for R packages
# devtools        | Tools to make package development easier
# lintr           | Static code analysis for R (linting)
# pkgdown         | Generate documentation websites for R packages
# roxygen2        | In-source documentation for R packages
# testthat        | Unit testing framework for R
# tryCatchLog     | Enhanced error logging for tryCatch workflows

# 3. Reporting & Publishing
# blogdown        | Create blogs and websites with R Markdown
# bookdown        | Author books and technical documents with R Markdown
# kableExtra      | Enhanced table styling for R Markdown outputs
# knitcitations   | Citation management for R Markdown
# knitr           | Dynamic report generation engine
# magick          | Advanced image processing tools
# rmarkdown       | Create dynamic documents, reports, and presentations
# tinytex         | Lightweight LaTeX distribution for PDF rendering

# 4. Data Wrangling & Analysis
# car             | Companion to Applied Regression tools
# data.table      | Fast data manipulation and aggregation
# dplyr           | Grammar of data manipulation
# FSA             | Fisheries Stock Analysis functions
# gmodels         | Tools for model fitting and analysis
# janitor         | Data cleaning helpers
# lme4            | Linear and generalized mixed-effects models
# lubridate       | Date and time parsing/manipulation
# MASS            | Functions and datasets from Modern Applied Statistics with S
# MESS            | Miscellaneous statistical tools
# multcomp        | Multiple comparison procedures
# nlme            | Linear and nonlinear mixed-effects models
# readxl          | Read Excel files
# reshape2        | Flexible data reshaping
# rfishbase       | Interface to FishBase data
# tidymodels      | Modeling and machine learning framework
# tidyverse       | Core collection of data science packages
# vegan           | Community ecology analysis
# xtable          | Export tables to LaTeX/HTML
# plyr            | Data splitting, applying, and combining
# abind           | Combine multidimensional arrays by binding along dimensions
# reshape         | Data reshaping utilities (legacy reshape framework)

# 5. Visualization
# cowplot         | Plot themes and annotations for ggplot2
# fmsb            | Radar charts and miscellaneous plotting functions
# gganimate       | Animation for ggplot2 graphics
# ggbeeswarm      | Beeswarm-style categorical scatter plots
# ggforce         | Extensions for ggplot2
# ggmap           | Spatial visualization with ggplot2
# ggpattern       | Pattern fills for ggplot2 geoms
# ggpubr          | Publication-ready ggplot2 wrappers
# ggtext          | Improved text rendering in ggplot2
# ggthemes        | Extra themes and scales for ggplot2
# patchwork       | Combine multiple ggplots into one layout
# plotly          | Interactive web graphics
# plotrix         | Additional plotting functions
# png             | Read and write PNG images
# ragg            | High-quality graphics device
# rgl             | 3D visualization
# RColorBrewer    | Color palettes for maps and figures
# wesanderson     | Wes Anderson-inspired color palettes
# scales          | Scale functions and formatters for visualization

# 6. Special/Experimental
# htmlwidgets     | Interactive JavaScript visualizations from R
# shiny           | Web application framework for R
# shinyBS         | Bootstrap components for Shiny
# shinyjs         | JavaScript helpers for Shiny apps
# XLConnect       | Read/write/format Excel files (Java-based)
# xlsx            | Read/write Excel files (Java-based)

# 7. Alaska Work
# sf              | Simple Features support for vector spatial data
# terra           | Spatial raster/vector data processing
# raster          | Raster data manipulation and analysis
# easyNCDF        | Read and write NetCDF data
# DirichletReg    | Dirichlet regression models for compositional data

# 8. Project-Specific / GitHub Packages
# mizer             | Size spectrum ecological modeling
# mizerExperimental | Experimental features for mizer
# mizerMR           | Multi-species and management-rule extensions for mizer workflows
# mizerReef         | Reef ecosystem extensions for mizer
# akgfmaps          | Mapping utilities for Alaska fisheries workflows
# coldpool          | Cold pool and temperature index calculations
# akfishcondition   | Fish condition metrics for Alaska species
# catchfunction     | Catch functions for size spectrum and fisheries modeling

# System Maintenance
# installr        | Update R and Rtools (Windows only)