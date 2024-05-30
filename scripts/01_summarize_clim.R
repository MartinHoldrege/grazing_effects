# Purpose: use monthly climate data from STEPWAT2 to calculate seasonality
# metrics

# Author: Martin Holdrege

# Started: May 14, 2024


# dependencies ------------------------------------------------------------

library(tidyverse)
source('src/paths.R')
source('src/general_functions.R')
# read in data ------------------------------------------------------------
# file KP provided on 5/12/2024. 
#This file contains monthly precip and temp for each of 150 years for all 200 sites for all current and future climate scenarios
clim1 <- read_csv(file.path(path_sw, 'Grazing.Project.2024.temp.precip.short.csv'),
                  col_types = 'iiiddicccciiicc')


# prepare dataframe -------------------------------------------------------


clim1 <- clim1 %>% 
  dplyr::select(Year, Month, TEMP_avg_C_Mean, PRECIP_ppt_Mean, site, GCM, years, RCP)


# summarize  ------------------------------------------------

# mean proportion summer precip (mean of proportions, 
# same way as calculated for Holdrege et al. 2024, Fire Ecology)
psp <- clim1 %>% 
  group_by(site, GCM, years, RCP, Year) %>% 
  summarize(psp = sum(PRECIP_ppt_Mean[Month %in% c(6, 7, 8)])/
              sum(PRECIP_ppt_Mean),
            # correlation between monthly t and p, calculated fo each year
            PTcor2 = stats::cor(TEMP_avg_C_Mean, PRECIP_ppt_Mean),
            AP = sum(PRECIP_ppt_Mean),
            AT = mean(TEMP_avg_C_Mean),
            .groups = 'drop_last') %>% 
  summarize(psp = mean(psp),
            PTcor2 = mean(PTcor2),
            MAP = mean(AP),
            MAT = mean(AT))


avg1 <- clim1 %>% 
  group_by(site, GCM, years, RCP, Month) %>% 
  summarize(temp = mean(TEMP_avg_C_Mean),
            ppt = mean(PRECIP_ppt_Mean), 
            .groups = 'drop_last')
# precipitation seasonality (correlationg between mean monthly ppt and temp)
seas <- avg1 %>% 
  summarize(PTcor = cor(temp, ppt))

out <- seas %>% 
  left_join(psp) %>% 
  mutate(RCP = rcp2factor(RCP),
         years = years2factor(years))


# save output -------------------------------------------------------------

write_csv(out, 'data_processed/site_means/clim_seasonality.csv')
