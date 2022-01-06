# Martin Holdrege

# Script started: 1/4/2022

# Purpose: create useful summary dataframes based on the.
# csv(s) created in the 01_query_db.R script
# this script will then be sourced by later scripts to create figures etc.


# dependencies ------------------------------------------------------------

library(tidyverse)
source("functions/general_functions.R")

# read in files -----------------------------------------------------------

# site level means of biomass across years, for each treatment and PFT
# created in 01_query_db.R script
bio3 <- read_csv("data_processed/site_means/bio_mean_by_site-PFT.csv",
                 show_col_types = FALSE)


# parse -------------------------------------------------------------------
# Convert columns to useful factors

bio4 <- bio3 %>% 
  mutate(graze = graze2factor(intensity),
         years = years2factor(years),
         RCP = rcp2factor(RCP)) %>% 
  dplyr::select(-intensity)


# note: when 0 biomass is simulated the PFT still shows up, but just with 
# 0 biomass

# total biomass by PFT ----------------------------------------------------
# creating dataframes where biomass is grouped into PFT's of interest

group_cols <- c('years', 'RCP', 'graze', 'site', "PFT", 'GCM')
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
# STOP: Continue here--diff should be done for each GCM then median
# calculated
pft5_bio_d1 <-  pft5_bio1 %>% 
  group_by(site, PFT, graze) %>%  
  scaled_change(by = c("PFT", "graze")) %>% 
  # median across GCMs
  group_by(site, years, RCP, PFT, graze) %>% 
  summarise(biomass = median(biomass),
            bio_diff = median(bio_diff),
            .groups = "drop") 



