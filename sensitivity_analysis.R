#title: "Sensitivity Analysis"
#format: html
#editor: visual

# Clear workspace and free up memory
# rm(list = ls())
# gc()

# Load required packages
library(tidyverse);
library(ggplot2);
library(readxl);
library(lubridate);
library(dplyr);
library(scales);

here::i_am("sensitivity_analysis.R")

# Load sensitivity data
final_scores <- read_excel(here("data/Sensitivity_Final_Scores.xlsx"), sheet = "Sensitivity",
                           col_names = c("scorer", "workshop", "species", "attribute_name", "low", "moderate", "high", "very_high", "tallies", "data_quality"),
                           skip = 1, na = "")

final_scores <- final_scores |> mutate(
  group = workshop |> replace_values(c("Rockfish1", "Rockfish2") ~ "Rockfish",
                                    c("Flatfish1", "Flatfish2") ~ "Flatfish",
                                    "Forage" ~ "Forage Fish")
)

## replace all NAs with 0
final_scores[is.na(final_scores)] <- 0

final_scores$attribute_name <- gsub("Sensitivity to Ocean Acidification", "Sensitivity to OA", final_scores$attribute_name)
final_scores$species <- gsub("Rougheye/blackspotted rockfish", "Rougheye_blackspotted rockfish", final_scores$species)

## recode attributes 
final_scores <- final_scores |>  mutate(
  attribute_name = factor(attribute_name),
  scorer = factor(scorer),
  species = factor(species),
  group = factor(group),
  workshop = factor(workshop),
  attribute_short = recode(attribute_name, 
                           "Habitat Specificity" = "HS",
                           "Thermal Tolerance" = "TT", 
                           "Sensitivity to OA" = "OA", 
                           "Foraging Strategy" = "FS", 
                           "Adult Movement" = "AM", 
                           "Dispersal Capability" = "DC", 
                           "Parental Investment" = "PI", 
                           "Reproductive Plasticity" = "RP", 
                           "Spawning Duration" = "SD", 
                           "Life History Strategy" = "LH", 
                           "Stock Status" = "SS", 
                           "Genetic Diversity" = "GD"),
  attribute_short = factor(attribute_short)
) 

final_scores <- final_scores |>
  subset(tallies > 0)
  

final_scores$attribute_name = as.factor(final_scores$attribute_name)
final_scores$attribute_name = ordered(final_scores$attribute_name, 
                                      levels = c("Habitat Specificity",
                                                 "Thermal Tolerance",
                                                 "Sensitivity to OA",
                                                 "Foraging Strategy",
                                                 "Adult Movement",
                                                 "Dispersal Capability", 
                                                 "Parental Investment",
                                                 "Reproductive Plasticity", 
                                                 "Spawning Duration",
                                                 "Life History Strategy",
                                                 "Stock Status", 
                                                 "Genetic Diversity"))


# ============================================================================
# SECTION 1: Calculate Sensitivity Scores
# ============================================================================
# Sums tallies across scorers, calculates weighted means for each sensitivity
# attribute, applies logic model (Morrison et al. 2015), calculates sensitivity 
# scores by group, and makes relevant plots.
# ============================================================================

# function 1
calculate_sum_scores <- function(data) {
  # sum tallies for each species
  calculated_sum_scores <- data |>
    group_by(workshop, group, species, attribute_short, attribute_name) |>
    summarize(sum_low = sum(low), sum_moderate = sum(moderate), sum_high = sum(high), sum_vh = sum(very_high)
    ) |>
    mutate( # calculate average score for each sensitivity attribute by species
      weighted_mean = ((sum_low*1)+(sum_moderate*2)+(sum_high*3)+(sum_vh*4))/(sum_low+sum_moderate+sum_high+sum_vh),
    )
  
  return(calculated_sum_scores)
}

# function 2
calculate_sensitivity_scores <- function(data) {
  sum_scores <- calculate_sum_scores(data)
  # apply logic model
  calculated_sensitivity_scores <- sum_scores |>
    group_by(workshop, group, species) |>
    summarize(
      above_3.5 = sum(weighted_mean >= 3.5),
      above_3 = sum(weighted_mean >= 3),
      above_2.5 = sum(weighted_mean >=2.5),
      sensitivity = ifelse(above_3.5 >= 3, "very high", ifelse(above_3 >= 2, "high", ifelse(above_2.5 >= 2, "moderate", "low"))),
    ) |>
    ungroup() |>
    select(!c(above_3.5, above_3, above_2.5))
  
  calculated_sensitivity_scores$sensitivity = ordered(calculated_sensitivity_scores$sensitivity, 
                                           levels = c("low",
                                                      "moderate",
                                                      "high",
                                                      "very high"))
  
  return(calculated_sensitivity_scores)
}

alternative_sensitivity_scores <- function(data) {
  sum_scores <- calculate_sum_scores(data)
  # apply logic model
  calculated <- sum_scores |>
    group_by(workshop, group, species) |>
    summarize(
      above_3.5 = sum(weighted_mean >= 3.5),
      above_2.5 = sum(weighted_mean >= 2.5), # mean scores between 2.5 and 3.5 would be considered as “high"
      above_1.5 = sum(weighted_mean >= 1.5), # mean scores between 1.5 and 2.5 would be considered as “moderate”
      sensitivity = ifelse(above_3.5 >= 3, "very high", ifelse(above_2.5 >= 2, "high", ifelse(above_1.5 >= 2, "moderate", "low"))),
    ) |>
    ungroup() |>
    dplyr::select(!c(above_3.5, above_2.5, above_1.5))
  
  calculated$sensitivity = ordered(calculated$sensitivity, 
                                                      levels = c("low",
                                                                 "moderate",
                                                                 "high",
                                                                 "very high"))
  
  return(calculated)
}

# calculations from functions 1 and 2
# sum_scores <- calculate_sum_scores(final_scores)
sensitivity_scores <- calculate_sensitivity_scores(final_scores)
alternative_sensitivity_scores <- alternative_sensitivity_scores(final_scores)
sensitivity_scores <- sensitivity_scores |>
  mutate(sensitivity_number = recode(sensitivity, 
                                  "low" = "1",
                                  "moderate" = "2", 
                                  "high" = "3", 
                                  "very high" = "4"))
alternative_sensitivity_scores <- alternative_sensitivity_scores |>
  mutate(sensitivity_number = recode(sensitivity, 
                                     "low" = "1",
                                     "moderate" = "2", 
                                     "high" = "3", 
                                     "very high" = "4"))

# save as csv files for vulnerability analysis
write.csv(sensitivity_scores, "sensitivity_scores_original.csv", row.names = FALSE)
write.csv(alternative_sensitivity_scores, "sensitivity_scores_alternative.csv", row.names = FALSE)


group_scores <- sensitivity_scores |>
  summarize(count = n(), .by = c(sensitivity, group)) |>
  left_join(sensitivity_scores |> summarize(total = n(), .by = c(group))) |>
  mutate(group = paste(group, " (n = ", total, ")", sep = "")) |>
  mutate(prop = count / total) 

group_scores$sensitivity = ordered(group_scores$sensitivity,
                                      levels = c("low",
                                                 "moderate",
                                                 "high",
                                                 "very high"))

# calculate proportion of tallies
sum_scores <- sum_scores |>
  mutate(
    total_tallies = sum_low + sum_moderate + sum_high + sum_vh,
    prop_low = sum_low/total_tallies,
    prop_moderate = sum_moderate/total_tallies,
    prop_high = sum_high/total_tallies,
    prop_vh = sum_vh/total_tallies
  )

# distribution plots of sensitivity scoring for each species (for workshops)
workshop_plots <- function(data, species_name) {
  ## summary of tallies barplot
  filtered <- data |>
    filter(species == species_name) |>
    pivot_longer(cols=c("prop_low", "prop_moderate", "prop_high", "prop_vh"),
                 names_to = "category",
                 values_to = "score") |>
    mutate(
      category_short = recode(category, "prop_low" = "L", "prop_moderate" = "M", "prop_high" = "H", "prop_vh" = "VH"),
      category_short = factor(category_short, levels=c("L", "M", "H", "VH"))
    )
  
  plot <- ggplot(filtered) +
    geom_col(mapping = aes(x = category_short, y = score, fill = category_short), position = "stack", linewidth = 0.25, colour="black", show.legend = FALSE) +
    xlab("Sensitivity Attribute") +
    ylab("Percentage of Tallies") +
    scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0,1),
      breaks = seq(0, 1, by = 0.2),
    ) +
    facet_wrap(~attribute_name) +
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, size = 10),
          strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  setwd("~/Documents/OSU/GOA CVA/Sensitivity/Final Scores/Species Plots/")
  ggsave(paste(species_name, "sensitivity.png", sep="_"), plot = plot, width=8, height=5)
  
}

# for loop to create all species plots
species_name<-unique(sum_scores$species)
for (i in 1:length(species_name)) {
  workshop_plots(sum_scores, species_name[i])
}

# group comparison plots
setwd("~/Documents/OSU/GOA CVA/Sensitivity/Final Scores/Group Plots")
ggplot(group_scores) +
  geom_col(mapping = aes(x = sensitivity, y = prop, fill = sensitivity), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
  facet_wrap(vars(group)) +
  scale_x_discrete(labels = c("L", "M", "H", "V")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of Group") +
  xlab("Sensitivity Score") +
  scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
  theme_bw() + 
  theme(strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
ggsave(filename="group_sensitivity_scores.png", plot = get_last_plot(), device = "png", width = 7, height = 5, bg = "transparent", dpi = 300)

# ============================================================================
# SECTION 2: Calculate Potential for Distributional Shift
# ============================================================================
# Use a subset of sensitivity attributes to evaluate the likelihood that a species
# will exhibit a shift in distribution with climate change (high adult movement,
# high dispersal capability, low habitat specificity, and low thermal tolerance). 
# Follows same methods (weighted means, logic model, etc.) as sensitivity and 
# creates plot to compare distribution shift potential by functional group.
# ============================================================================

# distribution shift potential
# based on four attributes: high adult movement, high dispersal capability, 
# low habitat specificity, low thermal tolerance

distributional_shift <- sum_scores |> # 
  filter(attribute_short == "HS" | attribute_short == "DC" | attribute_short == "AM" | attribute_short == "TT") |> # filter HS, AM, and DC
  mutate(
    dist_low = ifelse(attribute_short == "HS" | attribute_short == "AM" | attribute_short == "DC", sum_vh, sum_low),
    dist_moderate = ifelse(attribute_short == "HS" | attribute_short == "AM" | attribute_short == "DC", sum_high, sum_moderate),
    dist_high = ifelse(attribute_short == "HS" | attribute_short == "AM" | attribute_short == "DC", sum_moderate, sum_high),
    dist_vh = ifelse(attribute_short == "HS" | attribute_short == "AM" | attribute_short == "DC", sum_low, sum_vh),
    dist_weighted_mean = ((dist_low*1)+(dist_moderate*2)+(dist_high*3)+(dist_vh*4))/(dist_low+dist_moderate+dist_high+dist_vh),
  ) |>
  select(!c(sum_low, sum_moderate, sum_high, sum_vh, weighted_mean, prop_low, prop_moderate, prop_high, prop_vh))

dist_sensitivity_scores <- distributional_shift |>   # apply logic model
    group_by(group, species) |>
    summarize(
      above_3.5 = sum(dist_weighted_mean >= 3.5),
      above_3 = sum(dist_weighted_mean >= 3),
      above_2.5 = sum(dist_weighted_mean >=2.5),
      dist_potential = ifelse(above_3.5 >= 3, "very high", ifelse(above_3 >= 2, "high", ifelse(above_2.5 >= 2, "moderate", "low"))),
      dist_number = ifelse(above_3.5 >= 3, 4, ifelse(above_3 >= 2, 3, ifelse(above_2.5 >= 2, 2, 1))),
    ) |>
    ungroup() |>
    select(!c(above_3.5, above_3, above_2.5))

dist_sensitivity_scores$dist_potential <- factor(dist_sensitivity_scores$dist_potential, levels = c("low", "moderate", "high", "very high"))
dist_sensitivity_scores$dist_potential = ordered(dist_sensitivity_scores$dist_potential,
                                                 levels = c("low",
                                                            "moderate",
                                                            "high",
                                                            "very high"))
  
dist_group_scores <- dist_sensitivity_scores |>
    summarize(count = n(), .by = c(dist_potential, group)) |>
    left_join(dist_sensitivity_scores |> summarize(total = n(), .by = c(group))) |>
    mutate(group = paste(group, " (n = ", total, ")", sep = "")) |>
    mutate(prop = count / total) |>
    complete(dist_potential, group, fill = list(count = 0, prop = 0))

dist_group_scores$dist_potential <- factor(dist_group_scores$dist_potential, levels = c("low", "moderate", "high", "very high"))
dist_group_scores$dist_potential = ordered(dist_group_scores$dist_potential,
                                   levels = c("low",
                                              "moderate",
                                              "high",
                                              "very high"))

# group comparison plots for distributional shift 
setwd("~/Documents/OSU/GOA CVA/Sensitivity/Final Scores/Distributional Shift")
ggplot(dist_group_scores) +
  geom_col(mapping = aes(x = dist_potential, y = prop, fill = dist_potential), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
  facet_wrap(vars(group)) +
  scale_x_discrete(labels = c("L", "M", "H", "V")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of Group") +
  xlab("Potential for Distribution Shift") +
  scale_fill_brewer(palette = "BuPu") +
  theme_bw() + 
  theme(strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
ggsave(filename="group_dist_shift.png", plot = get_last_plot(), device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)



# ============================================================================
# SECTION 3: Calculate Directional Effect Scores and Compare
# ============================================================================
# Import directional effect scores, calculate directional effect scores across 
# species and groups, and make comparison plots.
# ============================================================================

# Load directional effect data
de_scores <- read_excel(here("data/Sensitivity_Final_Scores.xlsx"), sheet = "Directional effect",
                           col_names = c("scorer", "workshop", "species", "negative", "neutral", "positive"),
                           skip = 1, na = "")

de_scores <- de_scores |> mutate(
  group = workshop |> replace_values(c("Rockfish1", "Rockfish2") ~ "Rockfish",
                                     c("Flatfish1", "Flatfish2") ~ "Flatfish",
                                     "Forage" ~ "Forage Fish")
)

de_scores[is.na(de_scores)] <- 0
de_scores$species <- gsub("Rougheye/blackspotted rockfish", "Rougheye_blackspotted rockfish", de_scores$species)
de_scores <- de_scores |>  mutate(
  scorer = factor(scorer),
  species = factor(species),
  group = factor(group))

# calculate average directional effect by species
de_scores <- de_scores |>  
  mutate(
    score = ((negative*-1)+(neutral*0)+(positive*1))/(negative+neutral+positive)
  )

average_de <- de_scores |>
  group_by(group, species) |>
  summarize(de_mean = mean(score))

average_de <- average_de |> mutate(
  directional_effect = ifelse(de_mean <= -0.333, "negative", ifelse(de_mean >= 0.333, "positive", "neutral"))
) |>
  ungroup()

de_group_scores <- average_de |>
  summarize(count = n(), .by = c(directional_effect, group)) |>
  left_join(average_de |> summarize(total = n(), .by = c(group))) |>
  mutate(prop = count / total) |>
  mutate(group = paste(group, " (n = ", total, ")", sep = "")) |>
  complete(directional_effect, group, fill = list(count = 0, prop = 0))

de_group_scores$directional_effect = ordered(de_group_scores$directional_effect,
                                           levels = c("negative",
                                                      "neutral",
                                                      "positive"))

# group comparison plots for directional effect
setwd("~/Documents/OSU/GOA CVA/Sensitivity/Final Scores/Directional Effect")
ggplot(de_group_scores) +
  geom_col(mapping = aes(x = directional_effect, y = prop, fill = directional_effect), position = "dodge", linewidth = 0.25, colour="black", width = 0.8, show.legend = FALSE) +
  facet_wrap(vars(group)) +
  scale_x_discrete(labels = c("Negative", "Neutral", "Positive")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Proportion of Group") +
  xlab("Directional Effect") +
  scale_fill_manual(values = c("red", "beige", "green")) +
  theme_bw() + 
  theme(strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
ggsave(filename="group_directional_effect.png", plot = get_last_plot(), device = "png",width = 7, height = 5, bg = "transparent", dpi = 300)

# join and export distribution and directional effect scores

additional_scores <- left_join(dist_sensitivity_scores, average_de |> select(-group), by = "species") 
write.csv(additional_scores, "additional_scores.csv", row.names = FALSE)

# ============================================================================
# SECTION 4: Leave-One-Out Analysis
# ============================================================================
# Re-run each sensitivity calculation leaving out each scorer. This will point to 
# individual scorers whose scores had a large effect on changing the overall 
# sensitivity score for each species.
# ============================================================================

# create data frame to add original and LOO scores into
loo_analysis <- final_scores |>
  group_by(scorer, species) |>
  summarize() |>
  left_join(sensitivity_scores, join_by(species)) |>
  rename(sensitivity_with = sensitivity) 
  
# loop over scorers:
#   get list of species that scorer scored
#   recalculate sensitivity scores w/o that scorer for those species
#   add those scores to a dataframe
for (i in 1:length(loo_analysis$scorer)){
  scorer_name <- loo_analysis[[i, "scorer"]]
  species_name <- loo_analysis[[i, "species"]]
  
  # calculate sensitivity score, filtering to remove a particular scorer's score and include only one species
  # this creates 1 row, and then you pick what is in the sensitivity column
  sensitivity_without <- calculate_sensitivity_scores(final_scores |> filter(scorer != scorer_name, species == species_name))[[1,"sensitivity"]]
  
  loo_analysis[i,"sensitivity_without"] <- sensitivity_without
}

## count number of participants in each workshop
num_participants <- final_scores |>
  group_by(scorer, workshop) |>
  summarize() |>
  ungroup() 
num_participants2 <- num_participants |> summarize(participants = n(), .by = c(workshop)) |> ungroup()
num_participants <- num_participants |> left_join(num_participants2, join_by(workshop))

# identify potential outliers
potential_outliers <- loo_analysis |>
  left_join(num_participants, join_by(scorer)) |>
  mutate(
    change = ifelse(sensitivity_with != sensitivity_without, TRUE, FALSE) # did the score change without the removed scorer?
  ) |> filter(change == TRUE) |>  # filter out unchanged scores
  filter(participants > 4) # filter for workshops with more than 4 participants

potential_outliers_unique <- potential_outliers |>
  group_by(species, sensitivity_with, sensitivity_without) |>
  mutate(dupe = n()>1) |>
  filter(dupe == FALSE)

# make box plot of average scores for each 

scorer_averages <- final_scores |>
  group_by(workshop, scorer, species, attribute_short, attribute_name) |>
  summarize(sum_low = sum(low), sum_moderate = sum(moderate), sum_high = sum(high), sum_vh = sum(very_high)
  ) |>
  mutate( # calculate average score for each sensitivity attribute by species and scorer
    weighted_mean = ((sum_low*1)+(sum_moderate*2)+(sum_high*3)+(sum_vh*4))/(sum_low+sum_moderate+sum_high+sum_vh),
  )

scorer_averages <- scorer_averages |>
  group_by(workshop, scorer, species) |>
  summarize(mean = mean(weighted_mean))

# function that makes a plot of the average participant score for each species in a given workshop
outlier_plot <- function(data, workshop_name) {
  filtered <- data |>
    filter(workshop == workshop_name)
  
  plot <- ggplot(filtered) +
    geom_boxplot(aes(x = species, y = mean)) +
    geom_point(aes(x = species, y = mean, color=scorer), size=2) +
    xlab("Species") +
    ylab("Average Participant Score") +
    scale_x_discrete(labels = wrap_format(10)) +
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, size = 10),
        strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  setwd("~/Documents/OSU/GOA CVA/Sensitivity/Final Scores/Scorer LOO/")
  ggsave(paste(workshop_name, "scorer_distribution.png", sep="_"), plot = plot, width=8, height=5)
}

outlier_plot(scorer_averages, "Groundfish")

# for loop to create boxplots for each workshop - don't use these for anything
workshop_name<-unique(scorer_averages$workshop)
for (i in 1:length(workshop_name)) {
  outlier_plot(scorer_averages, workshop_name[i])
}

# ============================================================================
# SECTION 5: Data quality scores
# ============================================================================
# Calculate average data quality score for each species.
# ============================================================================

data_quality_sensitivity <- final_scores |>
  summarize(mean = mean(data_quality), .by = c(species)) |>
  mutate(across(where(is.numeric), round, digits = 1))

# ============================================================================
# SECTION 6: Other plots
# ============================================================================
# Individual plots for each species' sensitivity attribute.
# ============================================================================

# # one attribute for one species plot
# ggplot(final_scores |> filter(species == "Dusky rockfish", attribute_short == "HS")) +
#   geom_col(mapping = aes(x = category_short, y = score, fill = category_short), linewidth = 0.25, colour="black", show.legend = FALSE) +
#   xlab("Scoring Category") +
#   ylab("Expert Score Distribution") +
#   scale_fill_manual(values = c("green", "yellow", "orange", "red")) +
#   scale_y_continuous(
#     limits = c(0,21),
#     breaks = seq(0, 20, by = 5), 
#   ) +
#   theme_bw() +
#   theme(strip.text = element_text(hjust = 0, size = 12),
#         strip.background = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
#         rect = element_rect(fill = "transparent", color = "transparent", linewidth = 0),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# ggsave("~/Documents/OSU/GOA CVA/Sensitivity/Final Scores/dusky_habitat.png", plot = last_plot(), width=8, height=5)

