# Martin Holdrege

# Script started: 1/4/2022

# Purpose: create useful summary dataframes based on the.
# csv(s) created in the 01_query_db.R script
# this script will then be sourced by later scripts to create figures etc.

# Notes:
# Consider calculating agreement metrics (proportion of
# GCMs that 'agree' on the direction of change for a given site)

# dependencies ------------------------------------------------------------

library(tidyverse)
source("src/general_functions.R")

# read in files -----------------------------------------------------------

# site level means of biomass across years, for each treatment and PFT
# created in 01_query_db.R script
bio3 <- read_csv("data_processed/site_means/bio_mean_by_site-PFT.csv",
                 show_col_types = FALSE)


# parse -------------------------------------------------------------------
# Convert columns to useful factors

bio4 <- bio3 %>% 
  mutate(graze = graze2factor(intensity),
         RCP = rcp2factor(RCP),
         years = years2factor(years),
         # the unique combination of treatments
         id = paste(RCP, years, graze, sep = "_"), 
         id = str_replace(id, " ", "")) %>% 
  dplyr::select(-intensity) %>% 
  # using arrange so that id is an appropriately ordered factor
  arrange(graze, RCP, years) %>% 
  mutate(id = factor(id, levels = unique(id)))


# note: when 0 biomass is simulated the PFT still shows up, but just with 
# 0 biomass

# total biomass by PFT ----------------------------------------------------
# creating dataframes where biomass is grouped into PFT's of interest

# GCM--needs to be listed last for sequential grouping below
group_cols <- c('years', 'RCP', 'graze', 'site', "PFT", "id", 'GCM')
# * pft5 ------------------------------------------------------------------

# 5 main pft categories
pft5_bio1 <- bio4 %>% 
  mutate(PFT = pft5_factor(PFT)) %>%
  group_by(across(all_of(group_cols))) %>% 
  # b/some PFTs combined
  summarize(biomass = sum(biomass), .groups = "drop_last") %>% 
  dplyr::filter(!is.na(PFT))

pft5_bio2 <- pft5_bio1 %>% 
  # median across GCMs
  summarise(biomass = median(biomass),
            .groups = "drop")

# % change in biomass by PFT ----------------------------------------------

# * pft5 ------------------------------------------------------------------

# d stands for 'difference'
# % change in biomass from current conditions ,
# scaled by maximum biomass under current conditions(for a given grazing trmt)

pft5_bio_d1 <-  pft5_bio1 %>% 
  group_by(site, PFT, graze) %>%  
  scaled_change(by = c("PFT", "graze")) %>% 
  # median across GCMs
  group_by(site, years, RCP, PFT, graze, id) %>% 
  summarise(biomass = median(biomass),
            bio_diff = median(bio_diff),
            .groups = "drop") 

#adding another id variable
pft5_bio_d2 <- pft5_bio_d1 %>% 
  arrange(RCP, years) %>% # for creating ordered factor
  # id variable grazing removed
  mutate(id2 = str_replace(id, "_[A-z]+$", ""),
         id2 = factor(id2, levels = unique(id2)))

# wildfire ----------------------------------------------------------------

fire1 <- bio4 %>% 
  # fire return interval. WildFire is the mean number of fires in a given year
  # across 200 iterations
  mutate(fire_return = 1/(WildFire/200)) %>% 
  # taking average for each plot (otherwise value is repeated for each PFT)
  group_by(across(all_of(group_cols[group_cols != "PFT"]))) %>% 
  summarize(fire_return = mean(fire_return), .groups = "drop_last") %>% 
  # median across GCMs
  summarize(fire_return = median(fire_return, na.rm = TRUE), .groups = "drop") %>% 
  # if no fires occurred then set fire return interval to NA
  mutate(fire_return = ifelse(is.infinite(fire_return), NA, fire_return))

# change in fire return interval
fire_d1 <- fire1 %>% 
  group_by(graze) %>% 
  # warning here is ok
  scaled_change(var = "fire_return", by = "graze")


