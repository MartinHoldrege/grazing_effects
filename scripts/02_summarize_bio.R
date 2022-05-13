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
library(dtplyr) # to speed up a couple slow chunks
source("src/general_functions.R")

# read in files -----------------------------------------------------------

# site level means of biomass across years, for each treatment and PFT
# created in 01_query_db.R script
bio3 <- read_csv("data_processed/site_means/bio_mean_by_site-PFT.csv",
                 show_col_types = FALSE)

# parse -------------------------------------------------------------------
# Convert columns to useful factors

bio4a <- bio3 %>% 
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

# sites where no c4 is simulated under current conditions (or any conditions
# for the c4off simulation)
sites_noc4 <- bio4a %>% 
  filter(c4 == "c4off", biomass == 0, PFT == "p.warm.grass", RCP == "Current") %>%
  pull(site) %>% 
  unique()

stopifnot(length(sites_noc4) == 98)


# * deal with indivs column -----------------------------------------------
# num of individuals is just given for the 4 letter species codes,
# not the written out functional types they correspond to

# GCM, followed by site--needs to be listed last for sequential grouping below
group_cols <- c("c4", 'years', 'RCP', 'graze',  "id", "PFT", 'site', 'GCM')

indivs1 <- bio4a %>% 
  select(-biomass) %>% 
  filter(!is.na(indivs)) %>% 
  # convert spp codes to PFT
  mutate(PFT = spp2pft(PFT)) %>% 
  select(all_of(group_cols), indivs)

# add indivs back into the main df
bio4 <- bio4a %>% 
  select(-indivs) %>% 
  left_join(indivs1, by = group_cols)


# climate -----------------------------------------------------------------
# current climate

clim1 <- bio4 %>% 
  # arbitrarily filtering for one graze and PFT, so rows aren't duplicated
  filter(years == "Current", graze == "Light", PFT == "sagebrush",
         c4 == "c4on") %>% 
  select(site, PPT, Temp)


# total biomass by PFT ----------------------------------------------------
# creating dataframes where biomass is grouped into PFT's of interest


# * pft5 ------------------------------------------------------------------
# for historical code reasons this is called pft5 (the main 5 pfts), but
# now this includes all PFT groupings, including pft5

# functions that label PFTs into various groupings
factor_funs <- list(pft5_factor, pft3_factor, Pgrass_factor, pft_total_factor)

pft5_bio0 <- map(factor_funs, function(f) {
  out <- bio4 %>% 
    # create factor
    mutate(PFT = f(PFT)) %>%
    group_by(across(all_of(group_cols))) %>% 
    # b/some PFTs combined
    summarize(biomass = sum(biomass), 
              indivs = sum(indivs),
              .groups = "drop_last") %>% 
    dplyr::filter(!is.na(PFT))
  out
})


# examine PFT levels created
map(pft5_bio0, function(df) unique(df$PFT))

# includes total biomass as a separate PFT
pft5_bio1_tot <- bind_rows(pft5_bio0) %>% 
  mutate(PFT = pft_all_factor(PFT))

# total biomass not included (currently the default for subsequent figures)
pft5_bio1 <- bind_rows(pft5_bio0) %>% 
  filter(PFT!= "Total") %>% 
  # order PFTs into a factor
  mutate(PFT = pft_all_factor(PFT))

pft5_bio2 <- pft5_bio1 %>% 
  # median across GCMs
  summarise(biomass = median(biomass),
            indivs = median(indivs),
            .groups = "drop")%>% 
  left_join(clim1, by = "site") # adding current climate

# for now we're exluding succulents--is that warranted?--it appears so,
# they have very limited biomass
x <- bio4 %>% filter(PFT == "succulents") %>% pull(biomass)
mean(x > 0) # proportion sites w/ some succulents
mean(x[x>0]) # mean biomass at sites where present
max(x)
mean(x) # mean succulent biomass

# ** c4 on vs off ---------------------------------------------------------

# wide format, of biomass and number of individuals
# w/ c4 on vs off only for sites were c4 not 
# simulated under current conditions
pft5_c4on_v_off <- pft5_bio2 %>% 
  filter(site %in% sites_noc4) %>% 
  pivot_wider_c4()

# check all 0's for c4 grass
stopifnot(with(pft5_c4on_v_off,
               c4off[PFT == "C4Pgrass"] == 0))

pft5_c4on_v_off <- pft5_c4on_v_off %>% 
  filter(PFT != "C4Pgrass")

# % change of going from c4on to c4off, this needs to be calculated
# at the GCM level first b/ it is not a linear transformation.
# This df includes total biomass
c4on_v_off_diff <- pft5_bio1_tot %>% 
  pivot_wider_c4() %>% 
  mutate(bio_diff = (c4off -c4on)/c4on*100, # % change
         bio_es = log(c4off/c4on),# effect size
         bio_es = ifelse(is.finite(bio_es), bio_es, NA),
         indivs_diff = (indivs_c4off -indivs_c4on)/indivs_c4on*100, # % change
         indivs_es = log(indivs_c4off/indivs_c4on),# effect size
         indivs_es = ifelse(is.finite(indivs_es), indivs_es, NA),) %>% 
  group_by(across(all_of(group_cols[!group_cols %in% c("c4", "GCM")]))) %>% 
  lazy_dt() %>%  # for speed
  summarise(across(.cols = c(bio_diff, bio_es, indivs_diff, indivs_es),
                   .fns = median, na.rm = TRUE),
            .groups = "drop") %>% 
  as_tibble()

# C3/Pgrass ---------------------------------------------------------------

# proportion of Pgrass biomass that is C3Pgrass
C3_Pgrass_ratio <- pft5_bio1 %>% 
  filter(PFT %in% c("C3Pgrass", "Pgrass")) %>% 
  select(-indivs) %>% # pivot doesn't work with other data column present
  pivot_wider(names_from = "PFT",
              values_from = "biomass") %>%  
  mutate(C3_Pgrass_ratio = C3Pgrass/Pgrass) %>% 
  # median across GCMs
  group_by(across(all_of(group_cols[group_cols != "PFT"]))) %>% 
  summarise(C3_Pgrass_ratio = median(C3_Pgrass_ratio),
            .groups = 'drop')

# % change in biomass by PFT ----------------------------------------------


# ** change relative to same grazing trmt -----------------------------------

# d stands for 'difference'
# % change in biomass and num individuals from current conditions ,
# scaled by maximum biomass under current conditions(for a given grazing trmt)


pft5_bio_d2 <-  pft5_bio1 %>% 
  # calculating % scaled change for biomass and individuals
  scaled_change_2var(by = c("c4", "PFT", "graze")) %>% 
  # median across GCMs
  group_by(c4, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_diff") %>% 
  create_id2() # adding another id variable


# effect size; es = ln(biomass scenario of interest/biomass references group)
# effect size compared to current scenario of the given grazing scenario

pft5_bio_es1 <- pft5_bio1 %>% 
  # warning here is ok
  scaled_change_2var(by = c("c4", "PFT", "graze"), percent = FALSE,
                effect_size = TRUE) %>% 
  # median across GCMs
  group_by(c4, site, years, RCP, PFT, graze, id) %>% 
  lazy_dt() %>% # so next line run by dtpyr (for speed)
  summarise_bio_indivs(suffix = "_es") %>% 
  create_id2() %>% 
  mutate(bio_es = ifelse(is.finite(bio_es), bio_es, NA_real_),
         indivs_es = ifelse(is.finite(indivs_es), indivs_es, NA_real_))


# ** change relative to reference graze ------------------------------------
# comparing all grazing treatments and RCPs/time periods to a given reference 
# grazing scenario (e.g. current light grazing). 
# here using pft5_bio2 (already median across GCMs), because in all cases
# comparing to current conditions so calculating % change or effect size
# for each GCM then taking the median would yield the same answer

levs_graze <- levels(pft5_bio1$graze)
names(levs_graze) <- levs_graze

# scaled % changes
# naming here: d== difference, grefs = different grazing references used
pft5_d_grefs <- map(levs_graze, function(x) {
  out <-  pft5_bio2 %>% # using data already summarized across GCMs
    scaled_change_2var(by = c("c4", "PFT"), ref_graze = x) %>% 
    # adding id variable that doesn't include graze
    create_id2()
  out
})

# effect sizes
pft5_es_grefs <- map(levs_graze, function(x) {
  out <-  pft5_bio2 %>% # using data already summarized across GCMs
    scaled_change_2var(by = c("c4", "PFT"), ref_graze = x, percent = FALSE,
                  effect_size = TRUE) %>% 
    create_id2()
  out
})


# ** change relative to light grazing of same gcm -------------------------


# e.g. this shows the effects size of going from light grazing, to heavy
# grazing for RCP 8.5 end of century

# naming" d = difference, wgcm = within gcm comparison
pft5_d_wgcm <- pft5_bio1 %>% 
  scaled_change_2var(by = c("c4", "PFT", "RCP", "GCM", "years"), 
                ref_graze = "Light", percent = TRUE, effect_size = FALSE,
                within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(c4, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_diff") %>% 
  create_id2()


# naming es = effects size wgcm = within gcm comparison
pft5_es_wgcm <- pft5_bio1 %>% 
  scaled_change_2var(by = c("c4", "PFT", "RCP", "GCM", "years"), 
                ref_graze = "Light", percent = FALSE, effect_size = TRUE,
                within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(c4, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_es") %>% 
  create_id2()

# ** change relative to heavy grazing of same gcm -------------------------

# naming" d = difference, wgcm = within gcm comparison
pft5_d_wgcm_heavy <- pft5_bio1 %>% 
  scaled_change_2var(by = c("c4", "PFT", "RCP", "GCM", "years"), 
                     ref_graze = "Heavy", percent = TRUE, effect_size = FALSE,
                     within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(c4, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_diff") %>% 
  create_id2()

# threshold ---------------------------------------------------------------


# * % below threshold --------------------------------------------------
# percent of sites that are above a specified biomass threshold under future 
# climate/grazing

pcent <- 0.05 # 5th percentile is the reference level

# reference lower limit threshold
# --we want to know how many sites are going below this 
# threshold
# 5th percentile under current light grazing for each PFT, at sites
# with >0 biomass
ref_threshold <- pft5_bio1 %>% 
  filter(RCP == "Current", graze == "Light",
         biomass > 0) %>% 
  group_by(c4, PFT) %>% 
  # lower limit of biomass, under today's conditions
  summarise(threshold = quantile(biomass, probs = pcent),
            .groups = "drop")

threshold1 <- pft5_bio1 %>% 
  left_join(ref_threshold,  by = c("c4", "PFT")) %>% 
  group_by(across(all_of(group_cols[group_cols != "site"]))) %>% 
  mutate(above= biomass > threshold) %>% 
  # % sites above threshold
  summarise(pcent_above = mean(above)*100,
            .groups = "drop_last") %>% 
  # median across GCMs
  summarise(pcent_above = median(pcent_above),
            .groups = "drop")

# * crossing threshold ----------------------------------------------------
# What is mildest grazing treatment that causes biomass to go below threshold. 
# That is for a given scenario (most importantly, current and RCP8.5 mid century),
# what is the mildest grazing scenario that causes a site's biomass to go
# below the threshold (or is it already below threshold)

thresh_min_graze1 <- pft5_bio2 %>% 
  left_join(ref_threshold,  by = c("c4", "PFT")) %>% 
  mutate(above= biomass > threshold) %>% 
  group_by(c4, years, RCP, site, PFT) %>% 
  # min_graze = 1 (means it was already below the threshold w/ light grazing
  # min_graze == 5 means it was above the threshold for all grazing trmts)
  summarise(min_graze = cross_threshold(graze = graze, above = above),
            .groups = "drop")



# Composition -------------------------------------------------------------
# calculating the proportion of biomass made up by the 5 most important
# plant functional types


comp1 <- pft5_bio1_tot %>% 
  filter(c4 == 'c4on') %>% 
#  filter_rcp_c4() %>% 
  # median across GCMs
  summarise(biomass = median(biomass),
            .groups = 'drop') %>% 
  group_by(c4, years, RCP, graze,  id, PFT) %>% 
  # mean biomass across sites
  summarise(biomass = mean(biomass),
            .groups = 'drop_last') %>% 
  # percent of total biomass
  mutate(bio_perc = biomass/biomass[PFT == "Total"] * 100) %>% 
  filter(PFT %in% c(pft5_factor(NULL, return_levels = TRUE), "Total")) 

# calculating biomass and % of total biomass in the 'other' category,
# that is not biomass from the main 5 PFTs
# all climate scenarios
comp2_all <- comp1 %>% 
  mutate(other_biomass = biomass[PFT == "Total"] - sum(biomass[PFT != "Total"]),
         other_perc = 100 - sum(bio_perc[PFT != 'Total']),
         biomass = other_biomass,
         bio_perc = other_perc,
         PFT = "Other") %>%
  select(-other_biomass, -other_perc) %>% 
  distinct() %>% 
  bind_rows(comp1) %>% 
  filter(PFT != "Total") %>% 
  mutate(PFT = factor(PFT,  
                      levels = c("Other", "Cheatgrass", "Pforb",  
                                 "C4Pgrass", 
                                 "C3Pgrass", "Sagebrush")))%>% 
  arrange(RCP, graze) 

# just current and RCP8.5 mid-century
comp2 <- comp2_all %>%  
  filter_rcp_c4()%>% 
  mutate(id = factor(id, levels = unique(id)))


# composition of just key herbaceous PFTs
comp_herb1 <- comp2 %>% 
  filter(PFT %in% c('C3Pgrass', 'C4Pgrass', 'Pforb', "Cheatgrass")) %>% 
  # this df is already grouped
  mutate(bio_perc = biomass/(sum(biomass))*100)

# test that the percents do sum to 100
test <- summarise(comp_herb1, test = sum(bio_perc))$test
stopifnot(all.equal(test, rep(100, length(test))))

# wildfire ----------------------------------------------------------------

# * return interval -------------------------------------------------------

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


# * change in interval ----------------------------------------------------

# Calculated change as absolute difference (not scaled % change),
# due to extreme max values
fire_d1 <- fire1 %>% 
  # warning here is ok, calculating the actual (absolute) change, not % change
  scaled_change(var = "fire_return", by = c("c4", "graze"),
                percent = FALSE) %>% 
  # adding id variable that doesn't include graze
  create_id2()


# misc --------------------------------------------------------------------

# examining the difference between sagebrush and total shrub categories:

# x1 <- pft5_bio1 %>% 
#   filter(PFT == "Sagebrush")
# 
# x2 <- pft5_bio1 %>% 
#   filter(PFT == "Shrub")
# names(x2)
# 
# inner_join(x1, x2, by = c("c4", "years", "RCP", "graze","site", "id", "GCM")) %>% 
#   mutate(diff = biomass.y - biomass.x) %>% 
#   pull(diff) %>% 
#   hist()
