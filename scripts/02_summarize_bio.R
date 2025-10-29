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
source("src/mapping_functions.R") # for calc_low/high functions

# params ------------------------------------------------------------------
source('src/params.R')
n_years <- 50 # number of years of data from simulations we're using (i.e. 
# years 101-150)

n_iter <- 200 # number of iterations run

# runs with all grazing levels
runs_graze <- c('NoC4Exp' = "fire1_eind1_c4grass0_co20_2503", 
                'default' = "fire1_eind1_c4grass1_co20_2503")
ref_graze = opt$ref_graze
# read in files -----------------------------------------------------------

# site level means of biomass across years, for each treatment and PFT
# created in 01_query_db.R script
bio3 <- read_csv("data_processed/site_means/bio_mean_by_site-PFT_v6.csv",
                 show_col_types = FALSE)

# climate seasonality (file created in 01_summarize_clim.R)
seas1 <- read_csv('data_processed/site_means/clim_seasonality.csv',
                  show_col_types = FALSE)

# file created by 00_query_weather_db.R
clim_db <- read_csv("data_processed/site_means/dbWeather_200sites.csv",
                    show_col_types = FALSE)

# parse -------------------------------------------------------------------
# Convert columns to useful factors

bio4a <- bio3 %>% 
  # the 2402 runs were done with the same parameters as 2311 (but moderate to vheavy graze), so renaming
  # so they're grouped together
  mutate(run = str_replace(run, '2402', '2311')) %>% 
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

test <- bio4a %>% 
  group_by(run, PFT, RCP, years, GCM, graze) %>% 
  summarize(n = n())

stopifnot(test$n == 200) # should 1 observation per site per scenario

# note: when 0 biomass is simulated the PFT still shows up, but just with 
# 0 biomass

# sites where no c4 is simulated under current conditions (or any conditions
# for the c4off simulation)
# sites_noc4 <- bio4a %>% 
#   filter(c4 == "c4off", biomass == 0, PFT == "p.warm.grass", RCP == "Current") %>%
#   pull(site) %>% 
#   unique()

# stopifnot(length(sites_noc4) == 98)


# GCM, followed by site--needs to be listed last for sequential grouping below
group_cols <- c("run", 'years', 'RCP', 'graze',  "id", "PFT", 'site', 'GCM')

bio4 <- bio4a 

# climate -----------------------------------------------------------------

seas2 <- seas1 %>% 
  select(-PTcor) 
clim_all1 <- bio4 %>% 
  # arbitrarily filtering for one graze and PFT, so rows aren't duplicated
  filter(graze == "Light", PFT == "sagebrush",
         run == runs_graze['default']) %>% 
  select(site, PPT, Temp, years, RCP, GCM) %>% 
  left_join(seas2, by = join_by(site, years, RCP, GCM))

# current climate only
clim1 <- clim_all1 %>% 
  filter(years == "Current") %>% 
  left_join(clim_db, by = c("site" = "Site_id")) %>% 
  select(site, PPT, Temp, CorrTP2, psp) 

clim_all2 <- clim_all1 %>% 
  select(-MAP, -MAT) %>%  # these just based on monthly averages so, averaging
# temp could be wrong
  rename(MAP = PPT,
         MAT = Temp)

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
              utilization = sum(utilization),
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
  mutate(PFT = pft_all_factor(PFT)) %>% 
  calc_aherb() %>% 
  # for later summarizing
  group_by(run, years, RCP, graze, id, site, PFT)

# select runs that have for which complete set of grazing levels exist
# (used for filter in this and other scripts)
runs_graze0 <- pft5_bio1 %>% 
  ungroup() %>% 
  select(run, graze) %>% 
  distinct() %>% 
  group_by(run) %>% 
  summarize(n = n()) %>% 
  filter(n == 4) %>% 
  pull(run)

# make sure selecting runs for which all grazing levs exist
stopifnot(runs_graze %in% runs_graze0) 

pft5_bio2 <- pft5_bio1 %>% 
  # median across GCMs
  summarise(biomass_low = calc_low(biomass),
            biomass_high = calc_high(biomass),
            biomass = median(biomass),
            indivs_low = calc_low(indivs),
            indivs_high = calc_high(indivs),
            indivs = median(indivs),
            utilization_low = calc_low(utilization),
            utilization_high = calc_high(utilization),
            utilization = median(utilization),
            .groups = "drop")%>% 
  left_join(clim1, by = "site") # adding current climate

# for now we're exluding succulents--is that warranted?--it appears so,
# they have very limited biomass
x <- bio4 %>% filter(PFT == "succulents") %>% pull(biomass)
mean(x > 0) # proportion sites w/ some succulents
mean(x[x>0]) # mean biomass at sites where present

# ** c4 on vs off ---------------------------------------------------------
# see pre July, 2024 commits for code chunk

# C3/Pgrass ---------------------------------------------------------------

# proportion of Pgrass biomass that is C3Pgrass
# commenting out because not currently using this
C3_Pgrass_ratio <- pft5_bio1 %>% 
  filter(PFT %in% c("C3Pgrass", "Pgrass")) %>% 
  select(-indivs, -utilization) %>% # pivot doesn't work with other data column present
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

# this throws a warning (which is ok) because utilization is 0 for some groups
pft5_bio_d2 <-  pft5_bio1 %>% 
  filter(run %in% runs_graze) %>% 
  # calculating % scaled change for biomass and individuals
  scaled_change_2var(by = c("run", "PFT", "graze"),
                     vars = c("biomass", "indivs", "utilization"),
                     # regular % change
                     percent = TRUE,
                     divide_by_max = FALSE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_diff", include_low_high = TRUE) %>% 
  create_id2() # adding another id variable
  


# effect size; es = ln(biomass scenario of interest/biomass references group)
# effect size compared to current scenario of the given grazing scenario

# shouldn't throw warnings
pft5_bio_es1 <- pft5_bio1 %>% 
  filter(run %in% runs_graze) %>% 
  # warning here is ok
  scaled_change_2var(by = c("run", "PFT", "graze"), percent = FALSE,
                     vars = c("biomass", "indivs", "utilization"),
                effect_size = TRUE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, PFT, graze, id) %>% 
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

# commenting out because currently we only have 'light' grazing anyway
pft5_d_grefs <- map(levs_graze, function(x) {
  out <-  pft5_bio2 %>% # using data already summarized across GCMs
    filter(run %in% runs_graze) %>% 
    scaled_change_2var(by = c("run", "PFT"), ref_graze = x) %>%
    # adding id variable that doesn't include graze
    create_id2()
  out
})

# effect sizes
pft5_es_grefs <- map(levs_graze, function(x) {
  out <-  pft5_bio2 %>% # using data already summarized across GCMs
    filter(run %in% runs_graze) %>% 
    scaled_change_2var(by = c("run", "PFT"), ref_graze = x, percent = FALSE,
                       vars = c("biomass", "indivs", "utilization"),
                       effect_size = TRUE) %>%
    create_id2()
  out
})


# ** change relative to reference grazing of same gcm -------------------------

# e.g. this shows the effects size of going from moderate grazing, to heavy
# grazing for RCP 8.5 end of century

# naming" d = difference, wgcm = within gcm comparison
# throws warnings b/ 0 utilization in some sites (ok)
pft5_d_wgcm <- pft5_bio1 %>% 
  scaled_change_2var(by = c("run", "PFT", "RCP", "GCM", "years"), 
                     vars = c("biomass", "indivs", "utilization"),
                ref_graze = "Light", percent = TRUE, effect_size = FALSE,
                within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_diff") %>% 
  create_id2()


# naming es = effects size wgcm = within gcm comparison
pft5_es_wgcm <- pft5_bio1 %>% 
  scaled_change_2var(by = c("run", "PFT", "RCP", "GCM", "years"), 
                     vars = c("biomass", "indivs", "utilization"),
                ref_graze = ref_graze, percent = FALSE, effect_size = TRUE,
                within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, PFT, graze, id) %>% 
  summarise_bio_indivs(suffix = "_es") %>% 
  create_id2()

# ** change relative to heavy grazing of same gcm -------------------------
# no other grazing intensities currently available
# naming" d = difference, wgcm = within gcm comparison
# pft5_d_wgcm_heavy <- pft5_bio1 %>% 
#   scaled_change_2var(by = c("run", "PFT", "RCP", "GCM", "years"), 
#                      ref_graze = "Heavy", percent = TRUE, effect_size = FALSE,
#                      within_GCM = TRUE) %>% 
#   # median across GCMs
#   group_by(run, site, years, RCP, PFT, graze, id) %>% 
#   summarise_bio_indivs(suffix = "_diff") %>% 
#   create_id2()

pft5_es_wgcm_heavy <- pft5_bio1 %>%
  filter(run %in% runs_graze) %>% 
  scaled_change_2var(by = c("run", "PFT", "RCP", "GCM", "years"),
                     ref_graze = "Heavy", percent = FALSE, effect_size = TRUE,
                     vars = c("biomass", "indivs", "utilization"),
                     within_GCM = TRUE) %>%
  # median across GCMs
  group_by(run, site, years, RCP, PFT, graze, id) %>%
  summarise_bio_indivs(suffix = "_es") %>%
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
         biomass > 0, run %in% runs_graze) %>%
  group_by(run, PFT) %>%
  # lower limit of biomass, under today's conditions
  summarise(threshold = quantile(biomass, probs = pcent),
            .groups = "drop")

threshold1 <- pft5_bio1 %>%
  filter(run %in% runs_graze) %>% 
  left_join(ref_threshold,  by = c("run", "PFT")) %>%
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

# thresh_min_graze1 <- pft5_bio2 %>% 
#   left_join(ref_threshold,  by = c("run", "PFT")) %>% 
#   mutate(above= biomass > threshold) %>% 
#   group_by(c4, years, RCP, site, PFT) %>% 
#   # min_graze = 1 (means it was already below the threshold w/ light grazing
#   # min_graze == 5 means it was above the threshold for all grazing trmts)
#   summarise(min_graze = cross_threshold(graze = graze, above = above),
#             .groups = "drop")



# Composition -------------------------------------------------------------
# calculating the proportion of biomass made up by the 5 most important
# plant functional types
# not updated yet

comp1 <- pft5_bio1_tot %>%
  filter(run == runs_graze['default']) %>%
  # median across GCMs
  summarise(biomass = median(biomass),
            .groups = 'drop') %>%
  group_by(run, years, RCP, graze,  id, PFT) %>%
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
  filter_rcp_run()%>%
  mutate(id = factor(id, levels = unique(id)))


# composition of just key herbaceous PFTs
comp_herb1 <- comp2 %>%
  filter(PFT %in% c('C3Pgrass', 'C4Pgrass', 'Pforb', "Cheatgrass")) %>%
  # this df is already grouped
  mutate(bio_perc = biomass/(sum(biomass))*100)

# test that the percents do sum to 100
test <- summarise(comp_herb1, test = sum(bio_perc))$test
stopifnot(all.equal(test, rep(100, length(test))))

# total utilization ------------------------------------------------------

# total utilization across PFTs
util1 <- pft5_bio1 %>% 
  group_by(run, years, RCP, graze, id, site, GCM) %>% 
  summarize(utilization = sum(utilization),
            .groups = 'drop_last')

util_med1 <- util1 %>% 
  # median across GCMs
  summarise(utilization = median(utilization),
            .groups = "drop")%>% 
  left_join(clim1, by = "site") # adding current climate

# ** change relative to same grazing trmt -----------------------------------

# effect size; es = ln(utilization scenario of interest/util references group)
# effect size compared to current scenario of the given grazing scenario

# shouldn't throw warnings
util_es1 <- util1 %>% 
  filter(run %in% runs_graze) %>% 
  # warning here is ok
  scaled_change(by = c("run","graze"), percent = FALSE,
                     var = "utilization",
                     effect_size = TRUE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, 
           graze, id) %>% 
  summarise_bio_indivs(suffix = "_es",
                       col_names = 'utilization',
                       abbreviations = 'util') %>% 
  create_id2() 


# ** change relative to reference graze ------------------------------------
# comparing all grazing treatments and RCPs/time periods to a given reference 
# grazing scenario (e.g. current light grazing). 
# here starting with median across GCMs, because in all cases
# comparing to current conditions so calculating % change or effect size
# for each GCM then taking the median would yield the same answer

# scaled % changes
# naming here: d== difference, grefs = different grazing references used

# effect sizes
util_es_grefs <- map(levs_graze, function(x) {
  out <-  util_med1 %>% # using data already summarized across GCMs
    filter(run %in% runs_graze) %>% 
    scaled_change(by = c("run"), ref_graze = x, percent = FALSE,
                       var = 'utilization',
                       effect_size = TRUE) %>%
    create_id2()
  out
})


# ** change relative to light grazing of same gcm -------------------------

# e.g. this shows the effects size of going from light grazing, to heavy
# grazing for RCP 8.5 end of century


# naming es = effects size wgcm = within gcm comparison
util_es_wgcm <- util1 %>% 
  scaled_change(by = c("run", "RCP", "GCM", "years"), 
                     var = "utilization",
                     ref_graze = ref_graze, percent = FALSE, effect_size = TRUE,
                     within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, graze, id) %>% 
  summarise_bio_indivs(suffix = "_es",
                       col_names = 'utilization',
                       abbreviations = 'util') %>% 
  create_id2()

# ** change relative to heavy grazing of same gcm -------------------------

util_es_wgcm_heavy <- util1 %>%
  filter(run %in% runs_graze) %>% 
  scaled_change(by = c("run", "RCP", "GCM", "years"),
                     ref_graze = "Heavy", percent = FALSE, effect_size = TRUE,
                     var = "utilization",
                     within_GCM = TRUE) %>%
  # median across GCMs
  group_by(run, site, years, RCP, graze, id) %>%
  summarise_bio_indivs(suffix = "_es",
                       col_names = 'utilization',
                       abbreviations = 'util') %>%
  create_id2()

# wildfire ----------------------------------------------------------------

# * probability/return interval ------------------------------------------

# one row for each GCM
fire0 <- bio4 %>% 
  # fire return interval. WildFire is the mean number of fires in a given year
  # across 200 iterations
  filter(PFT == 'sagebrush', 
         run %in% runs_graze,
         str_detect(run, 'fire1') # only want simulations where fire was actually simulated
         ) %>% 
  select(-PFT) %>% 
  # number of fires over all iterations and years
  rename(n_fires = WildFire) %>% 
  # WildFire is the sum of wildfires across iterations and years
  mutate(fire_prob = n_fires/(n_iter*n_years)*100,
         fire_return = (n_years*n_iter)/n_fires,
         fire_return = ifelse(is.infinite(fire_return), NA, fire_return)) %>%  # annual wildfire
  select(all_of(group_cols[group_cols != "PFT"]), matches("fire"))


# median across GCMs
fire_med1 <- fire0 %>% 
  group_by(run, years, RCP, graze, id, site) %>% 
  # recalculating fire return here
  summarize(fire_prob = median(fire_prob, na.rm = TRUE),
            n_fires = median(n_fires),
            .groups = "drop") %>% 
  # recalculate fire return because taking the median screws up the math 
  # (because 'Infinite' fire returns are removed with na.rm = TRUE)
  mutate(fire_return = (n_years*n_iter)/n_fires,
         fire_return = ifelse(is.infinite(fire_return), NA, fire_return))



# * change in interval ----------------------------------------------------

# Calculated change as absolute difference (not scaled % change),
# due to extreme max values
# change due to to climate within a grazing level
fire_d_wgraze <- fire_med1 %>%
  # warning here is ok, calculating the actual (absolute) change, not % change
  scaled_change_2var(vars = c("fire_return", "fire_prob"), by = c("run", "graze"),
                percent = FALSE) %>%
  # adding id variable that doesn't include graze
  create_id2()


# * graze effect within GCM -----------------------------------------------

# e.g. this shows the effects of going from light grazing, to heavy
# grazing for RCP 8.5 end of century

# naming" d = difference, wgcm = within gcm comparison

fire_d_wgcm <- fire0 %>% 
  scaled_change(by = c("run", "RCP", "GCM", "years"), 
                     var = "fire_prob",
                     ref_graze = "Light", percent = FALSE, effect_size = FALSE,
                     within_GCM = TRUE) %>% 
  # median across GCMs
  group_by(run, site, years, RCP, graze, id) %>% 
  summarise(fire_prob_diff = median(fire_prob_diff)) %>% 
  create_id2()

stopifnot(sum(is.na(fire_d_wgcm$fire_prob_diff)) == 0)

# save outputs ------------------------------------------------------------

out <- list(
  pft5_bio1 = pft5_bio1,
  pft5_bio2 = pft5_bio2,
  pft5_bio_d2 = pft5_bio_d2,
  pft5_bio_es1 = pft5_bio_es1,
  pft5_d_wgcm = pft5_d_wgcm,
  pft5_d_grefs = pft5_d_grefs,
  pft5_es_wgcm = pft5_es_wgcm,
  fire0 = fire0,
  fire_med1 = fire_med1,
  fire_d_wgraze = fire_d_wgraze,
  fire_d_wgcm = fire_d_wgcm,
  clim_all2 = clim_all2,
  clim1 = clim1,
  runs_graze = runs_graze,
  group_cols = group_cols
)

saveRDS(out, 'data_processed/site_means/summarize_bio.RDS')
