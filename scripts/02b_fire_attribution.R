# Purpose: Determine which predictor variables of the fire equation
# are most responsible the change in fire probability at a given site
# under a given climate scenario (and fixed grazing)
# or the change relative to moderate grazing under a fixed climate (i.e. 'gref')

# Author: Martin Holdrege

# Started: March 19, 2025


# params ------------------------------------------------------------------

# if change less than 1 fire per 500 years then 
# don't consider the driver of change
delta_fire_prob_cutoff <- 0.2 # percentage point change
ref_graze <- 'Moderate'
# dependencies ------------------------------------------------------------

library(tidyverse)
source("src/fire_functions.R")
source("src/qm_functions.R")
source("src/general_functions.R")
source("src/fig_functions.R")

# read in data ------------------------------------------------------------

# list of dataframes created in "scripts/02_summarize_bio.R" 
l <- readRDS('data_processed/site_means/summarize_bio.RDS') 

clim_all2 <- l$clim_all2

bio1 <- l$pft5_bio1

# quantiles used in STEPWAT2 for quantile mapping in fire equation
# created in scripts/bio_matching/02_quantile_matching_simple.R
# NOTE--as the quantiles used are updated in STEPWAT, update the file that is read in
quants <- read_csv("data_processed/qm/quantiles_for_qm_0.5match_fire0_eind1_c4grass1_co20.csv")


# create qm functions -----------------------------------------------------
# functions for quantile mapping

qm_Aherb <- qm_quant_factory(
  from = quants$biomass_stepwat[quants$PFT == 'afgAGB'], 
  to = quants$biomass_rap[quants$PFT == 'afgAGB'])

qm_Pherb <- qm_quant_factory(
  from = quants$biomass_stepwat[quants$PFT == 'pfgAGB'], 
  to = quants$biomass_rap[quants$PFT == 'pfgAGB'])


# prepare dataframes  --------------------------------------------------------

bio2 <- bio1 %>% 
  filter(str_detect(run, 'fire1')) %>% 
  #calc_aherb() %>% 
  filter(PFT %in% c('Aherb', 'Pherb')) %>% 
  select(run, years, RCP, graze, id, site, GCM, biomass, PFT) %>% 
  left_join(clim_all2, by = join_by(years, RCP, site, GCM))

bio3a <- bio2 %>% 
  pivot_wider(values_from = biomass, names_from = PFT) %>% 
  df_factor()

# calculate the quantile mapped biomass,
# so that predictions made with fire equation are appropriate
bio3 <- bio3a %>% 
  mutate(Aherb = qm_Aherb(Aherb),
         Pherb = qm_Pherb(Pherb))

bio_cur1 <- bio3 %>% 
  filter(years == 'Current')

bio_fut1 <- bio3 %>% 
  filter(years != 'Current')

pred_vars <- c("MAP", "MAT", "psp", "Aherb", "Pherb")

bio_cur2 <- bio_cur1 %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_cur',
               values_to = "value_cur")  %>% 
  ungroup() %>% 
  select(-years, -RCP, -GCM, -id) 

# columns labeled with the pred vars, are the 'current' values
bio_cur3 <- bio_cur2 %>% 
  full_join(bio_fut1, by = c('run', 'graze', 'site'),
             relationship = "many-to-many")

# change in predictor variable (within a grazing level)
delta_pred_var1 <- bio_fut1 %>% pivot_longer(cols = all_of(pred_vars),
                     names_to = 'pred_var_fut',
                     values_to = "value_fut") %>% 
  left_join(bio_cur2, by = c('run', 'graze', 'site', 
                             'pred_var_fut' = 'pred_var_cur')) %>% 
  mutate(delta_pred_value = value_fut - value_cur) %>% 
  rename(pred_var = pred_var_fut)


test <- bio_cur3 %>% 
  filter(site == 1,
         years == levels(years)[3],
         RCP == levels(RCP)[2],
         GCM == unique(GCM)[1],
         run == unique(run)[1],
         graze == "Light")

stopifnot(nrow(test) == length(pred_vars))

# all pred vars (columns) are equal to future values, 
# except the one named in pred_var_cur column, where the current (historical)
# value is instead used
one_change <- bio_cur3 %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_fut',
               values_to = 'value_fut') %>% 
  # replacing one pred var at a time, with the current value
  mutate(value = ifelse(pred_var_fut == pred_var_cur, value_cur, value_fut)) %>% 
  select(-value_cur, -value_fut) %>% 
  pivot_wider(values_from = value, names_from = pred_var_fut) 


# * grazing reference -----------------------------------------------------
# here 'ref' refers to the reference grazing level, and 'grz' refers
# to any other grazing level. 

bio_ref1 <- bio3 %>% 
  filter(graze == ref_graze)

bio_grz <- bio3 %>% 
  filter(graze != ref_graze)


bio_ref2 <- bio_ref1 %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_ref',
               values_to = "value_ref")  %>% 
  ungroup() %>% 
  select(-graze,  -id) 

# columns labeled with the pred vars, are the 'reference' values
bio_ref3 <- bio_ref2 %>% 
  full_join(bio_grz, by = c('run', 'site', 'years', 'GCM', 'RCP'),
            relationship = "many-to-many")

# change in predictor variable (within a climate scenario)
delta_pred_var_gref1 <- bio_grz %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_grz',
               values_to = "value_grz") %>% 
  left_join(bio_ref2, by = c('run', 'site', 'years', 'GCM', 'RCP',
                             'pred_var_grz' = 'pred_var_ref')) %>% 
  mutate(delta_pred_value = value_grz - value_ref) %>% 
  rename(pred_var = pred_var_grz)


# all pred vars (columns) are equal to value for given grazing level, 
# except the one named in pred_var_ref column, where the
# value for the reference grazing level is used is instead used
one_change_gref <- bio_ref3 %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_grz',
               values_to = 'value_grz') %>% 
  # replacing one pred var at a time, with the current value
  mutate(value = ifelse(pred_var_grz == pred_var_ref, value_ref, value_grz)) %>% 
  select(-value_ref, -value_grz) %>% 
  pivot_wider(values_from = value, names_from = pred_var_grz) %>% 
  # only biomass variables are different between grazing levels 
  # (climate is fixed)
  filter(pred_var_ref %in% c('Aherb', 'Pherb'))



# predict fire probability ------------------------------------------------

predict_fire_with <- function(df) {
  p <- with(df, predict_fire(mat = MAT, map = MAP, psp = psp, afg = Aherb,
                        pfg = Pherb,run_checks = TRUE))
  p*100 # convert to %
}

# predicted fire probability when 1 variable changed to the current value
one_change$fire_prob_1var <- predict_fire_with(one_change)

# predicted fire probability under future conditions
bio_fut2 <- bio_fut1
bio_fut2$fire_prob_fut <- predict_fire_with(bio_fut2)

# predicted fire probability under current conditions
bio_cur1b <- bio_cur1
bio_cur1b$fire_prob_cur <- predict_fire_with(bio_cur1b)

# combine
one_change2 <- bio_cur1b %>% 
  ungroup() %>% 
  select(run, site, graze, fire_prob_cur) %>% 
  right_join(one_change) %>% 
  select(-all_of(pred_vars)) %>% 
  left_join(select(bio_fut2, run, years, RCP, GCM, graze, id, site, fire_prob_fut)) %>% 
  # How much different is fire probability, when 1 variable isn't changed
  # (i.e. 1 variable changed to current value)
  mutate(delta_1var = fire_prob_fut - fire_prob_1var ,
         rcp_years = rcp_label(RCP, years, include_parenth = FALSE)) 

one_change2b <- delta_pred_var1 %>% 
  rename(pred_value_fut = value_fut, pred_value_cur = value_cur,
         pred_var_cur = pred_var) %>% 
  right_join(one_change2, by = c("run", "years", "graze", "site",  "RCP", 
                                 "GCM", 'pred_var_cur', "id")) 


# * grazing reference ----------------------------------------------------

# predicted fire probability when 1 variable changed to the current value
one_change_gref$fire_prob_1var <- predict_fire_with(one_change_gref)

# predicted fire probability under future conditions
bio_grz2 <- bio_grz
bio_grz2$fire_prob_grz <- predict_fire_with(bio_grz2)

# predicted fire probability under current conditions
bio_ref1b <- bio_ref1
bio_ref1b$fire_prob_ref <- predict_fire_with(bio_ref1b)

# combine
one_change_gref2 <- bio_ref1b %>% 
  ungroup() %>% 
  select(run, site, RCP, GCM, years, fire_prob_ref) %>% 
  right_join(one_change_gref) %>% 
  select(-all_of(pred_vars)) %>% 
  left_join(select(bio_grz2, run, years, RCP, GCM, graze, id, site, fire_prob_grz)) %>% 
  # How much different is fire probability, when 1 variable isn't changed
  # (i.e. 1 variable changed to current value)
  mutate(delta_1var = fire_prob_grz - fire_prob_1var ,
         rcp_years = rcp_label(RCP, years, include_parenth = FALSE)) 

one_change_gref2b <- delta_pred_var_gref1 %>% 
  rename(pred_value_grz = value_grz, pred_value_ref = value_ref,
         pred_var_ref = pred_var) %>% 
  right_join(one_change_gref2, by = c("run", "years", "graze", "site",  "RCP", 
                                 "GCM", 'pred_var_ref', "id")) 


# summarize results -------------------------------------------------------

dominant_driver0 <- one_change2b %>% 
  rename(dominant_driver = pred_var_cur) %>% 
  # if change less than 1 fire per 500 years then 
  # don't consider the driver of change
  mutate(delta_fire_prob = fire_prob_fut - fire_prob_cur,
         dominant_driver = ifelse(abs(delta_fire_prob) < delta_fire_prob_cutoff, 
                                  'None', dominant_driver)) %>% 
  group_by(run, site, graze, RCP, years, id, GCM) %>% 
  filter(abs(delta_1var) == max(abs(delta_1var))) %>% 
  mutate(n = n()) 


# if n>1 that means there were ties for the dominant driver
# and these would need to be dealt with
stopifnot(all(dominant_driver0$n == 1))

# summary across GCMs
set.seed(1234)
dominant_driver1 <- dominant_driver0 %>% 
  group_by(run, site, graze, RCP, years, id) %>% 
  summarize(dominant_driver = mode_rand(dominant_driver),
            .groups = 'drop') %>% 
  mutate(dominant_driver = factor(dominant_driver,
                                  levels = c(!!pred_vars, 'None')))

# summarizie across GCMs
one_change3 <- one_change2b %>% 
  group_by(run, site, graze, years, RCP, id, pred_var_cur) %>% 
  summarize(across(.cols = c('fire_prob_fut', "fire_prob_cur", "delta_1var", 
                             "pred_value_fut", "pred_value_cur", "delta_pred_value"),
                   .fns = median),
            delta_1var_abs = median(abs(delta_1var)),
            .groups = 'drop') %>% 
  mutate(pred_var_cur = factor(pred_var_cur, levels = !!pred_vars))


# *grazing reference ------------------------------------------------------

dominant_driver_gref0 <- one_change_gref2b %>% 
  rename(dominant_driver = pred_var_ref) %>% 
  # if change less than 1 fire per 500 years then 
  # don't consider the driver of change
  mutate(delta_fire_prob = fire_prob_grz - fire_prob_ref,
         dominant_driver = ifelse(abs(delta_fire_prob) < delta_fire_prob_cutoff, 
                                  'None', dominant_driver)) %>% 
  group_by(run, site, graze, RCP, years, id, GCM) %>% 
  filter(abs(delta_1var) == max(abs(delta_1var))) %>% 
  mutate(n = n()) 


# if n>1 that means there were ties for the dominant driver
# and these would need to be dealt with
# but it is ok if the 'multiple' dominant drivers are all 'None'
test <- dominant_driver_gref0$dominant_driver[dominant_driver_gref0$n > 1]
stopifnot(all(test == 'None'))

# summary across GCMs
set.seed(1234)
dominant_driver_gref1 <- dominant_driver_gref0 %>% 
  group_by(run, site, graze, RCP, years, id) %>% 
  summarize(dominant_driver = mode_rand(dominant_driver),
            .groups = 'drop') %>% 
  mutate(dominant_driver = factor(dominant_driver,
                                  levels = c(!!pred_vars, 'None')))

tmp <- bio3a %>% 
  pivot_wider(id_cols = c('run', 'years', 'RCP', 'GCM', 'site'),
              values_from = all_of(pred_vars),
              names_from = 'graze')
with(tmp, plot(Aherb_Heavy, Aherb_Light))

# summarize across GCMs
one_change_gref3 <- one_change_gref2b %>% 
  group_by(run, site, graze, years, RCP, id, pred_var_ref) %>% 
  summarize(across(.cols = c('fire_prob_grz', "fire_prob_ref", "delta_1var", 
                             "pred_value_grz", "pred_value_ref", "delta_pred_value"),
                   .fns = median),
            delta_1var_abs = median(abs(delta_1var)),
            .groups = 'drop') %>% 
  mutate(pred_var_ref = factor(pred_var_ref, levels = !!pred_vars))

# checks ------------------------------------------------------------------

test <- dominant_driver1 %>% 
  group_by(run, id) %>% 
  summarize(n = n()) %>% 
  pull(n)
stopifnot(all(test == 200))  

test <- dominant_driver_gref1 %>% 
  group_by(run, id) %>% 
  summarize(n = n()) %>% 
  pull(n)
stopifnot(all(test == 200))  


# output ------------------------------------------------------------------

out <- list(
  one_change_smry1 = one_change3,
  one_change_gcm1 = one_change2b,
  one_change_gref_smry1 = one_change_gref3,
  one_change_gref_gcm1 = one_change_gref2b,
  dominant_driver1 = dominant_driver1,
  dominant_driver_gref1 = dominant_driver_gref1,
  pred_vars = pred_vars,
  pred_vars2 = levels(dominant_driver1$dominant_driver)
)

saveRDS(out, 'data_processed/site_means/fire_dominant_drivers.RDS')
