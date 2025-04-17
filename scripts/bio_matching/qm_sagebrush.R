# Purpose: quantile mapping for sagebrush stepwat biomass to RCMAP cover
# the cover was summarized the same way as perennial and annual
# biomass was for Holdrege et al (2024, fire ecology) which was used
# for quantile mapping (see scripts/bio_matching/01_combine_data.R).
# It was summarized the same way for consistency

# Note here using the latest interpolated stepwat run data, 
# because the

# parameters --------------------------------------------------------------

version_interp <- 'v4'
run <- "fire1_eind1_c4grass1_co20_2503"
# grazing levels to use
graze_levels <- c("L" = "Light", "M" = "Moderate", "H" = "Heavy")

qual_cutoff <- 0.5 # matching quality cutoff used determine how good
# the interpolation needed to be


run2 <- paste0(
  run, "_", paste0(names(graze_levels), collapse = "")
)
# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source('src/general_functions.R')

# read in data ------------------------------------------------------------

# * cover -----------------------------------------------------------------

# file created by the cheatgrass_fire/scripts/01_compile_sagebrush_ann.js script
# (https://code.earthengine.google.com/?scriptPath=users%2Fmholdrege%2Fcheatgrass_fire%3A01_compile_sagebrush_ann.js)
r_sage1 <- rast("data_raw/RCMAPV6_sagebrush-cover_1986-2019_3yrAvg_1000m_sagebrush-biome-mask_v2.tif")

# * interpolated stepwat --------------------------------------------------

pattern <- paste0(run, ".*Sagebrush.*Current.*.tif")

paths <- list.files(
  file.path("data_processed/interpolated_rasters/biomass", version_interp),
  pattern = pattern,
  full.names = TRUE)

basename(paths)
r_sw1 <- rast(paths)

# for filter which data to use
r_qual <- rast(paste0("data_processed/interpolation_quality/matching_quality",
                 version_interp, ".tif"))

# for later reproducibility (i.e., know which sites were used, where)
# gives the site number, of the site that was matched to the given pixel
r_site <- rast(paste0("data_processed/interpolation_data/interp_locations_200sites_",
                     version_interp, ".tif"))



# extract data ------------------------------------------------------------

cell_nums <- r_sage1[[1]]
names(cell_nums) <- 'cell_num'
cell_nums[] <- 1:ncell(r_sage1[[1]])
cell_nums[is.na(r_sage1[[1]])] <- NA

# ** extract rcmap data --------------------------------------------------------

rcmap_df1 <- as.data.frame(r_sage1)
rcmap_df1$cell_num <- as.numeric(row.names(rcmap_df1))

rcmap_df2 <- rcmap_df1 %>% 
  pivot_longer(cols = matches("rangeland_sagebrush"),
               values_to = "sagebrush_cover") %>% 
  mutate(year = as.numeric(str_extract(name, "\\d{4}$"))) %>% 
  select(-name)
  

# * interpolated stepwat data -------------------------------

# identify stepwat grid-cells w/ good matching, and that are also cells
# included in the fire manuscript
cell_nums2 <- extend(cell_nums, r_qual)
r_qual2 <- extend(r_qual, cell_nums2)

mask0 <- r_qual2 < qual_cutoff & !is.na(cell_nums2)
mask <- classify(mask0, rcl = cbind(FALSE, NA))
cell_nums_m <- mask(cell_nums2, mask)

r_site2 <- extend(r_site, cell_nums2) %>% 
  mask(mask)
names(r_site2) <- 'site'
r_sw2 <- extend(r_sw1, cell_nums2) %>% 
  mask(mask) %>% 
  c(., r_site2, cell_nums_m)



cell_nums_mv <- as.numeric(cell_nums_m[!is.na(cell_nums_m)]) # vector

rcmap_df_m <- rcmap_df2 %>% 
  filter(cell_num %in% cell_nums_mv) 

df_sw1 <- as_tibble(r_sw2) 

info1 <- create_rast_info(r_sw1)

df_sw2 <- df_sw1 %>% 
  pivot_longer(cols = -c(cell_num, site),
               names_to = 'id',
               values_to = 'biomass') %>% 
  left_join(info1, by = 'id') %>% 
  mutate(fire = str_extract(run, 'fire\\d')) %>% 
  # remove masked cells
  filter(graze %in% graze_levels)


# calculate quantiles -----------------------------------------------------

# important, that data from rcmap and stepwat
# corresponds the same grid-cells
stopifnot(rcmap_df_m$cell_num %in% df_sw2$cell_num,
          df_sw2$cell_num %in% rcmap_df_m$cell_num)


# quantiles to use
# (using more than needed, so don't need to intelligently
# select a smaller number)
probs <- c(0, 0.001, .005, 
            seq(0.01, 0.99, by = 0.02),
            0.995, 0.999, 1)

# for plotting
deciles <- tibble(
  deciles = seq(0, 1, by = 0.1),
  biomass_stepwat = quantile(df_sw2$biomass, probs = deciles)
)

quants <- tibble(quantile = probs)
quants$biomass_stepwat <- quantile(df_sw2$biomass, probs = probs)
quants$cover_rcmap <- quantile(rcmap_df_m$sagebrush_cover, probs = probs)

quants <- quants %>% 
  mutate(biomass_stepwat = signif(as.numeric(biomass_stepwat), digits = 3),
         cover_rcmap = signif(as.numeric(cover_rcmap)), digits = 3)

png(paste0("figures/bio_matching/q-q_plots_sagebrush_qm", qual_cutoff, "_simple_",run2, ".png"),
    width = 5, height = 5, units = 'in', res = 600)
ggplot(quants, aes(biomass_stepwat, cover_rcmap)) +
  geom_line() +
  geom_rug(data = deciles, aes(x = biomass_stepwat, y = NULL)) +
  labs(subtitle = 'Sagebrush biomass and cover, q-q plot, \nthese values used for quantile mapping',
       caption = paste0('rug shows deciles\n', run2)) +
  theme_bw()
dev.off()



# save output -------------------------------------------------------------

name <- paste0('quantiles_for_qm_sagebrush_', qual_cutoff, 'match_', run2, '.csv')

write_csv(quants, file.path('data_processed/qm', name))



