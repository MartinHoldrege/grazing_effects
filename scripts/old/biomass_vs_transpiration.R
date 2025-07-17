# Purpose--compare biomass and transpiration


# Author: Martin Holdrege
# Started: 9/27/2024


# dependencies ------------------------------------------------------------

library(tidyverse)


# params ------------------------------------------------------------------

run <- "fire1_eind1_c4grass0_co20_2311"
RCP <- 'Current'
graze <- c('grazeM' = 'Moderate')

# load data ---------------------------------------------------------------

# created in scripts/examine_transpiration.R
suffix <- paste0(str_replace_all(run, "(_)|(2311)", ""), names(graze), RCP)
t_list1 <- readRDS(paste0(
  'data_processed/site_means/transpiration_dfs_',
  suffix))

# created in scripts/02_summarize_bio.R
bio1 <- readRDS('data_processed/site_means/pft5_bio2.RDS')

# filter data -------------------------------------------------------------

t_by_pft1 <- t_list1$t_by_pft 

bio2 <- bio1 %>% 
  filter(.data$run == !!run,
         .data$RCP == !!RCP,
         .data$graze == !!graze)
  
# sum to grass, shrub, forb, to compare to soilwat2 T data
# lookup
bio3 <- bio2 %>% 
  mutate(PFT_sw = case_when(
    PFT == "Sagebrush" ~ 'shrub',
    PFT == 'Shrub' ~ 'shrub',
    PFT == 'Pgrass' ~ 'grass',
    PFT == 'Cheatgrass' ~ 'grass',
    PFT == 'Pforb' ~ 'forbs',
    PFT == 'Aforb' ~ 'forbs',
    TRUE ~ 'other'
  )) %>% 
  filter(PFT_sw != 'other') %>% 
  group_by(PFT_sw, site) %>% 
  summarize(biomass = sum(biomass),
            .groups = 'drop') %>% 
  rename(PFT = PFT_sw)

bio4 <- bio3 %>% 
  left_join(t_by_pft, by = c('PFT', 'site'))


# figures -----------------------------------------------------------------

pdf(paste0('figures/water/transpiration_vs_bio_', suffix, ".pdf"),
    width = 7, height = 6)
ggplot(bio4, aes(transp, biomass, color = MAT > 10)) +
  geom_point() +
  facet_wrap(~PFT, scales = 'free', ncol = 2) +
  geom_smooth(method = 'lm',
              se = FALSE) +
  labs(y = lab_bio0,
       x = lab_transp0,
       caption = 'transpiration values are from SOILWAT2, 
       but those are re-partitioned in STEPPE')
dev.off()

