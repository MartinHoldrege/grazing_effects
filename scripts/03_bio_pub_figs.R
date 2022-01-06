# Martin Holdrege

# Script started 1/6/2022

# Purpose--to create figures of biomass data from simulations, that are 
# good quality (i.e. goal is for them to be 'publication quality').
# This script does not create maps of the interpolated data. 

# dependencies ------------------------------------------------------------

library(tidyverse)
library(lemon) # for facet_rep_wrap function
source("scripts/02_summarize_bio.R")
source("src/fig_params.R")


# params ------------------------------------------------------------------

theme_set(theme_classic())
theme_update(strip.background = element_blank())
line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot

outlier.size = 0.5
# boxplots ----------------------------------------------------------------

# * absolute biomass ------------------------------------------------------
# This figure meant to be analogous to M.E.'s thesis figure 9.
# biomass by pft, rcp, time period, and grazing intensity

# dataframe of annotation for grazing levels in fig
graze_anno1 <- pft5_bio1 %>% 
  group_by(PFT, graze) %>% 
  # adding small 2.5% so the text fits over the main body of the figure
  summarize(y = max(biomass) + (max(biomass) - min(biomass))*0.025,
            x = median(as.numeric(id)),
            .groups = "drop_last") %>% 
  mutate(y = max(y)) # so all text is at same level in a given panel

jpeg("figures/biomass/pub_qual/pft5_bio_boxplot.jpeg",
     res = 600, height = 8, width = 5, units = "in")

pft5_bio2 %>% 
  # so boxplot doesn't show sites with 0 biomass
  mutate(biomass = ifelse(biomass == 0, NA, biomass)) %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  # first plotting text, so it doesn't overplot data
  geom_text(data = graze_anno1, 
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

# ** pft5 ------------------------------------------------------------------

# boxplot of change in biomass (scaled percent), for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period

# dataframe for geom_text
yr_anno1 <- pft5_bio_d2 %>% 
  group_by(PFT, RCP, years) %>% 
  summarise(x = median(as.numeric(id2)),
            y = max(bio_diff) + (max(bio_diff) - min(bio_diff))*0.15) %>%
  group_by(PFT) %>% 
  mutate(y = max(y),
         years = str_replace(years, "-", "-\n"))

jpeg("figures/biomass/pub_qual/pft5_bio_diff_boxplot.jpeg",
     res = 600, height = 8, width = 5, units = "in")

ggplot(pft5_bio_d2, aes(id2, bio_diff, fill = graze)) +
  geom_text(data = yr_anno1, 
            aes(x, y, label = years, fill = NULL),
            size = 2.5, vjust = "top", lineheight = .7) +
  geom_boxplot(position = "dodge",
               outlier.size = outlier.size) +
  facet_rep_wrap(~PFT, scales = "free", ncol = 2) +
  scale_fill_manual(values = cols_graze, 
                    name = "Grazing Treatment") +
  # so just display the RCP
  scale_x_discrete(labels = function(x) str_extract(x, "RCP(4.5|8.5)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = c(0.75, 0.15)) +
  geom_vline(xintercept = 2.5, linetype = 2) +
  labs(y = lab_bio2,
       x = lab_rcp)

dev.off()
