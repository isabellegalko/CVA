# ============================================================================
# Apply Bias Correction to All SSP Scenarios
# ============================================================================
# Purpose: Apply delta correction to all three SSP scenarios (126, 245, 585)
# and create a combined time series for analysis
# ============================================================================

# Clear workspace
rm(list = ls())
gc()

# Load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, arrow, viridis)
here::i_am("create_bias_corrected_ts.R")

# Source the bias correction function
source(here("functions.R"))

# ============================================================================
# Load Data
# ============================================================================

# turned this into a function to run for several ROMS variables!

bias_correction <- function(depth, variable_name) {
# Load hindcast and historical (same for all scenarios)
hindcast <- open_dataset(here("data/processed/hindcast_annual_data.parquet")) %>%
  filter(layer == depth, variable == variable_name) %>%
  collect()

historical <- open_dataset(here("data/processed/historical_annual_data.parquet")) %>%
  filter(layer == depth, variable == variable_name) %>%
  collect()

# Load all three SSP scenarios
ssp126 <- open_dataset(here("data/processed/ssp126_annual_data.parquet")) %>%
  filter(layer == depth, variable == variable_name) %>%
  collect()

ssp245 <- open_dataset(here("data/processed/ssp245_annual_data.parquet")) %>%
  filter(layer == depth, variable == variable_name) %>%
  collect()

ssp585 <- open_dataset(here("data/processed/ssp585_annual_data.parquet")) %>%
  filter(layer == depth, variable == variable_name) %>%
  collect()

# ============================================================================
# Apply Bias Correction to All Scenarios
# ============================================================================

# SSP1-2.6
cat("  Processing SSP1-2.6...\n")
corrected_ssp126 <- bias_correct_roms(
  hindcast = hindcast,
  historical = historical,
  projection = ssp126,
  use_sd = FALSE,
  include_hindcast = TRUE
)

# SSP2-4.5
cat("  Processing SSP2-4.5...\n")
corrected_ssp245 <- bias_correct_roms(
  hindcast = hindcast,
  historical = historical,
  projection = ssp245,
  use_sd = FALSE,
  include_hindcast = TRUE
)

# SSP5-8.5
cat("  Processing SSP5-8.5...\n")
corrected_ssp585 <- bias_correct_roms(
  hindcast = hindcast,
  historical = historical,
  projection = ssp585,
  use_sd = FALSE,
  include_hindcast = TRUE
)

# ============================================================================
# Combine Time Series
# ============================================================================
# The historical (1980-Jan 1990) and hindcast (Feb 1990-2020) portions are
# identical across all scenarios. We only need the projection portions to differ.

# Keep historical/hindcast from ssp126 only - it's conserved across runs
projection_126 <- corrected_ssp126 %>%
  mutate(run = case_when(date < as.Date("1990-02-01") ~ "historical",
                         date < as.Date("2021-01-01") ~ "hindcast",
                         .default = "ssp126"))

projection_245 <- corrected_ssp245 %>%
  filter(year > 2020) %>%
  mutate(run = "ssp245")

projection_585 <- corrected_ssp585 %>%
  filter(year > 2020) %>%
  mutate(run = "ssp585")

# Combine all
all_corrected <- bind_rows(
  projection_126,
  projection_245,
  projection_585
) %>%
  arrange(date)

# Save as parquet for efficient storage and future use
write_parquet(all_corrected, paste("data/processed/all_scenarios_bias_corrected_", variable_name, "_", depth, ".parquet", sep = ""))

cat(paste("Dataset saved to: data/processed/all_scenarios_bias_corrected_", variable_name, "_", depth, ".parquet\n", sep = ""))

}

bias_correction("surface", "temp")
bias_correction("bottom", "temp")
bias_correction("surface", "salt")
bias_correction("bottom", "salt")
bias_correction("surface", "PhL") # large phytoplankton concentration
bias_correction("surface", "PhS") # small phytoplankton concentration

bias_correction("surface", "Cop") # small copepod concentration
bias_correction("surface", "NCa") # large copepod concentration
bias_correction("surface", "Eup") # euphausiid concentration
bias_correction("surface", "MZL") # large microzooplankton concentration
bias_correction("surface", "MZS") # small microzooplankton concentration

# ============================================================================
# Visualize: Domain-Wide January Time Series
# ============================================================================

# Calculate January mean for each year and scenario
january_ts <- all_corrected %>%
  filter(month == 1) %>%
  group_by(year, run) %>%
  summarise(
    mean_temp = mean(value_dc, na.rm = TRUE),
    sd_temp = sd(value_dc, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
p1 <- ggplot(january_ts, aes(x = year, y = mean_temp, color = run)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_vline(xintercept = 2020, linetype = "dashed", color = "gray30") +
  annotate("text", x = 2020, y = max(january_ts$mean_temp), 
           label = "Scenarios\ndiverge", hjust = -0.1, size = 3, color = "gray30") +
  scale_color_manual(
    values = c(
      "historical" = "#1b9e77",
      "hindcast" = "#764",
      "ssp126" = "#2166ac",
      "ssp245" = "#ffa500",
      "ssp585" = "#b2182b"
    ),
    labels = c(
      "historical" = "Historical",
      "hindcast" = "Hindcast",
      "ssp126" = "SSP1-2.6 (Low)",
      "ssp245" = "SSP2-4.5 (Moderate)",
      "ssp585" = "SSP5-8.5 (High)"
    )
  ) +
  labs(
    title = "Bias-Corrected Bottom Temperature Projections for Gulf of Alaska",
    x = "Year",
    y = "Mean January Temperature (°C)",
    color = "Scenario",
    caption = "Historical and hindcast spliced through 2020; projections diverge after 2020"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

print(p1)

# ============================================================================
# Visualize: Monthly values model-wide
# ============================================================================

# Calculate monthly mean for each scenario (all months)
monthly_ts <- all_corrected %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(year_month, run) %>%
  summarise(mean_temp = mean(value_dc, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(monthly_ts, aes(x = year_month, y = mean_temp, color = run)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_vline(xintercept = as.Date("2020-12-31"), linetype = "dashed", color = "gray30") +
  scale_color_manual(
    values = c(
      "historical" = "#1b9e77",
      "hindcast" = "#764",
      "ssp126" = "#2166ac",
      "ssp245" = "#ffa500",
      "ssp585" = "#b2182b"
    ),
    labels = c(
      "historical" = "Historical",
      "hindcast" = "Hindcast",
      "ssp126" = "SSP1-2.6 (Low)",
      "ssp245" = "SSP2-4.5 (Moderate)",
      "ssp585" = "SSP5-8.5 (High)"
    )
  ) +
  labs(
    title = "Bias-corrected Monthly Bottom Temperature: All Scenarios",
    x = "Date",
    y = "Mean Temperature (°C)",
    color = "Scenario"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p2)

