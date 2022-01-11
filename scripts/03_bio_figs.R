# Martin Holdrege

# Script started 1/6/2022

# Purpose--to create figures of biomass data from simulations, that are 
# good quality 
# This script does not create maps of the interpolated data. 

# Note--many of these boxplots have repeated code, some of that could be 
# consolidated into functions

# dependencies ------------------------------------------------------------

library(tidyverse)
library(lemon) # for facet_rep_wrap function
source("scripts/02_summarize_bio.R")
source("src/fig_params.R")
source("src/fig_functions.R") # box_ann function defined here

# params ------------------------------------------------------------------

theme_set(theme_classic())
theme_update(strip.background = element_blank())
line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot

outlier.size = 0.5

# boxplots ----------------------------------------------------------------

# * absolute biomass ------------------------------------------------------
# This figure meant to be analogous to M.E.'s thesis figure 9.
# biomass by pft, rcp, time period, and grazing intensity

jpeg("figures/biomass/pub_qual/pft5_bio_boxplot.jpeg",
     res = 600, height = 8, width = 5, units = "in")

pft5_bio2 %>% 
  # so boxplot doesn't show sites with 0 biomass
  mutate(biomass = ifelse(biomass == 0, NA, biomass)) %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  # first plotting text, so it doesn't overplot data
  geom_text(data = ~box_anno(., var = "biomass", group_by = c("PFT", "graze"),
                             mult = 0.05),
            aes(x, y, label = graze, fill = NULL),
            size = 2.5) +
  geom_boxplot(outlier.size = outlier.size) + # not showing outliers as points
  facet_rep_wrap(~ PFT, scales = "free", ncol = 2) + 
  scale_fill_manual(values = cols_rcp, name = "Scenario") +
  scale_x_discrete(labels = years2lab) +
  geom_vline(xintercept = line_loc, linetype = 2) +
  theme(legend.position = c(0.75, 0.15),
        axis.text = element_text(size = 7)) +
  labs(x = lab_yrs,
       y = lab_bio0)
dev.off()


# * biomass change -------------------------------------------------------

# ** change relative to same graze ------------------------------------

# boxplot of change in biomass (scaled percent), for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period

jpeg("figures/biomass/pub_qual/pft5_bio_diff_boxplot.jpeg",
     res = 600, height = 8, width = 5, units = "in")

ggplot(pft5_bio_d2, aes(id2, bio_diff, fill = graze)) +
  geom_text(data = ~box_anno(., var = "bio_diff", group_by = c("PFT", "RCP"),
                             id = "id2", mult = 0.05),
            aes(x, y, label = RCP, fill = NULL),
            size = 2.5) +
  geom_boxplot(position = "dodge",
               outlier.size = outlier.size) +
  facet_rep_wrap(~PFT, scales = "free", ncol = 2) +
  scale_fill_manual(values = cols_graze, 
                    name = "Grazing Treatment") +
  # so just display the RCP
  scale_x_discrete(labels = years2lab) +
  theme(legend.position = c(0.75, 0.15)) +
  geom_vline(xintercept = 2.5, linetype = 2) +
  labs(y = lab_bio2,
       x = lab_yrs)

dev.off()

# ** change relative to reference graze ------------------------------------
# boxplots showing change in biomass relative to current time period and
# given grazing intensity

pdf("figures/biomass/pft5_bio_diff_gref_boxplots.pdf", 
    height = 6.5, width = 5)

map2(pft5_d_grefs, names(pft5_d_grefs), function(df, ref_graze){

  # at this point it is not possible to show the unused grazing level in ggplot
  # for the current rcp (https://github.com/tidyverse/ggplot2/issues/3345)
  g <- ggplot(df, aes(id2, bio_diff, fill = graze)) +
    geom_text(data = ~box_anno(., var = "bio_diff", 
                               group_by = c("PFT", "RCP"),
                               id = "id2", mult = 0.05),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5) +
    geom_hline(yintercept = 0, alpha = 0.3, linetype = 1) +
    geom_boxplot(position = position_dodge(preserve = "single"),
                 outlier.size = outlier.size) +
    facet_rep_wrap(~PFT, scales = "free", ncol = 2) +
    scale_fill_manual(values = cols_graze, 
                      name = "Grazing Treatment") +
    # so just display the RCP
    scale_x_discrete(labels = years2lab) +
    theme(legend.position = c(0.75, 0.15)) +
    geom_vline(xintercept = c(1.5, 3.5), linetype = 2) +
    labs(y = lab_bio2,
         x = lab_yrs,
         subtitle = paste("Change in biomass relative to", tolower(ref_graze), 
                          "grazing \n under current conditions"))
g
})
dev.off()


# fire --------------------------------------------------------------------

# boxplot of fire return interval by RCP, and grazing intensity
g1 <- fire1 %>% 
  group_by(id) %>% 
  mutate(ylim = boxplot.stats(fire_return)$stats[5]) %>% 
  filter(fire_return <= ylim) %>% 
  ggplot(aes(x = id, y = fire_return, fill = RCP)) +
  # first plotting text, so it doesn't overplot data
  geom_text(data = ~box_anno(., var = "fire_return",
                             group_by = c("graze")),
            aes(x, y, label = graze, fill = NULL),
            size = 2.5) +
  geom_boxplot(outlier.color = NA) + # not showing outliers as points
  scale_fill_manual(values = cols_rcp, name = "Scenario") +
  scale_x_discrete(labels = years2lab) +
  geom_vline(xintercept = line_loc, linetype = 2) +
  theme(legend.position = "top",
        axis.text = element_text(size = 7))  +
  labs(x = lab_yrs,
       y = lab_fire0)

# fire difference boxplot
g2 <- fire_d1 %>% 
  group_by(id) %>% 
  # removing outliers (extreme outliers make the body of the boxplot
  # hard to see)
  mutate(ymax = boxplot.stats(fire_return_diff)$stats[5],
         ymin = boxplot.stats(fire_return_diff)$stats[1],
         fire_return_diff = ifelse(fire_return_diff <= ymax &
                                     fire_return_diff >= ymin,
                                   fire_return_diff, NA_real_)) %>%
  ggplot(aes(id2, fire_return_diff, fill = graze)) +
  geom_text(data = ~box_anno(., var = "fire_return_diff", 
                             group_by = c("RCP"), id = "id2"),
            aes(x, y, label = RCP, fill = NULL),
            size = 2.5) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_boxplot(position = "dodge",
               outlier.color = NA) +
  scale_fill_manual(values = cols_graze, 
                    name = "Grazing Treatment") +
  # so just display the RCP
  scale_x_discrete(labels = years2lab) +
  theme(legend.position = "top") +
  geom_vline(xintercept = 2.5, linetype = 2) +
  labs(x = lab_yrs,
       y = lab_fire1,
       caption = "Outliers not shown in either panel")

jpeg("figures/fire/fire_return_boxplots.jpeg",
     height = 8, width = 5, res = 600, units = "in")
gridExtra::grid.arrange(
  g1, g2
)
dev.off()
