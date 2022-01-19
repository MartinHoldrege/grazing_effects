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

# for looping
levs_pft <- levels(pft5_bio2$PFT)

levs_pft_c4 <- expand_grid(pft = factor(levs_pft, levs_pft), 
                           levs_c4 = unique(pft5_bio2$c4)) %>% 
  arrange(pft, desc(levs_c4))

levs_c4 <- unique(levs_pft_c4$levs_c4) %>% 
  sort(., decreasing = TRUE)

# colors
scale_color_graze <- function() {
  scale_color_manual(values = cols_graze, name = "Grazing")
}

scale_fill_graze <- function() {
  scale_fill_manual(values = cols_graze, name = "Grazing")
}

# functions ---------------------------------------------------------------

# given a ggplot (g) create scatter plots y ~ MAP and y ~ MAT
climate_scatter <- function(g) {
  out <- list()
  # vs MAP
  out[["MAP"]] <- g +
    geom_point(aes(x = PPT)) +
    geom_smooth(aes(x = PPT), method = "loess", se = FALSE) +
    labs(x = lab_map)
  
  # vs MAT
  out[["MAT"]] <- g +
    geom_point(aes(x = Temp)) +
    geom_smooth(aes(x = Temp), method = "loess", se = FALSE) +
    labs(x = lab_mat)
  out
}


# absolute biomass ------------------------------------------------------
# This figure meant to be analogous to M.E.'s thesis figure 9.
# biomass by pft, rcp, time period, and grazing intensity


# * boxplot ---------------------------------------------------------------


# base of the next two figures
box1 <- function() {
  # first plotting text, so it doesn't overplot data
  list(
    geom_text(data = ~box_anno(., var = "biomass", group_by = c("PFT", "graze"),
                           mult = 0.05),
          aes(x, y, label = graze, fill = NULL),
          size = 2.5),
    geom_boxplot(outlier.size = outlier.size), # not showing outliers as points
    facet_rep_wrap(~ PFT, scales = "free", ncol = 2), 
    scale_fill_manual(values = cols_rcp, name = "Scenario"),
    scale_x_discrete(labels = years2lab),
    geom_vline(xintercept = line_loc, linetype = 2),
    theme(legend.position = c(0.75, 0.15),
          axis.text = element_text(size = 7)),
    labs(x = lab_yrs,
         y = lab_bio0)
  )
  }

jpeg("figures/biomass/pub_qual/pft5_bio_boxplot_c4on.jpeg",
     res = 600, height = 8, width = 5, units = "in")

pft5_bio2 %>% 
  filter(c4 == "c4on") %>% 
  # so boxplot doesn't show sites with 0 biomass
  mutate(biomass = ifelse(biomass == 0, NA, biomass)) %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  box1()

dev.off()

jpeg("figures/biomass/pub_qual/pft5_bio_boxplot_c4off.jpeg",
     res = 600, height = 8, width = 5, units = "in")

pft5_bio2 %>% 
  filter(c4 == "c4off") %>% 
  # so boxplot doesn't show sites with 0 biomass
  mutate(biomass = ifelse(biomass == 0, NA, biomass)) %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  box1()

dev.off()

# * scatterplot (vs climate) ----------------------------------------------

# each combination PFt and C4 on or off is plotted on seperate pages
# facets are RCP/year combinations, colors are grazing
pdf("figures/biomass/bio_vs_climate_v1.pdf",
    width = 6, height = 5)

# here pmap is working by taking the column names (i.e. iterating over rows)
pmap(levs_pft_c4, function(pft, levs_c4) {

  g <- pft5_bio2 %>% 
    filter(PFT == pft, c4 == levs_c4) %>% 
    ggplot(aes(y = biomass, color = graze)) +

    facet_rep_wrap(~RCP + years) +
    scale_color_graze() +
    labs(y = lab_bio0,
         subtitle = paste(pft, "biomass"),
         caption = c4on_off_lab(levs_c4)) +
    theme(legend.position = c(0.85, 0.15),
          axis.text = element_text(size = 7)) 

  climate_scatter(g)
  
})

# facets are grazing level, colors are RCP
pmap(levs_pft_c4, function(pft, levs_c4) {
  g <- pft5_bio2 %>% 
    filter(PFT == pft, c4 == levs_c4,
           years != "2030-2060") %>% 
    ggplot(aes(y = biomass, color = RCP)) +
    facet_rep_wrap(~graze) +
    scale_color_manual(values = cols_rcp) +
    labs(y = lab_bio0,
         subtitle = paste(pft, "biomass"),
         caption = paste(c4on_off_lab(levs_c4),
           "\nOnly showing 2070-2100 for RCP 4.5 & 8.5")) +
    theme(legend.position = "top",
          axis.text = element_text(size = 7))
  
  climate_scatter(g)
  
})

dev.off()


# * c4 on vs off biomass ----------------------------------------------------

pdf("figures/biomass/bio_c4on_vs_off_v1.pdf",
    width = 9, height = 6)
map(levs_pft[levs_pft != "C4Pgrass"], function(x) {
  pft5_c4on_v_off %>% 
    filter(PFT == x) %>% 
    ggplot(aes(c4on, c4off, color = graze)) +
    geom_point(alpha = 0.6, size = 1) +
    facet_rep_wrap(RCP ~ years) +
    geom_abline(slope = 1, color = "black") +
    geom_smooth(method = "lm") +
    labs(
      x = lab_c4on0,
      y = lab_c4off0,
      caption = paste(
        "Regression lines seperate for each grazing treatment,",
        "1:1 line in black.", 
      "Only showing the 98 sites for which no C4Pgrass is simulated\n",
      "under future conditions when the C4Pgrass climate suitability function", 
      "is turned off"),
      subtitle = paste(x, "biomass at the site level, for simulations with",
                       "C4Pgrass on vs off")
    ) +
    scale_color_graze() +
    theme(legend.position = c(0.85, 0.15))
})
dev.off()

# biomass change -------------------------------------------------------

# * change relative to same graze ------------------------------------

# boxplot of change in biomass (scaled percent), for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period

box2 <- function(){
  list(
    geom_text(data = ~box_anno(., var = "bio_diff", group_by = c("PFT", "RCP"),
                               id = "id2", mult = 0.05),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5),
    geom_boxplot(position = "dodge",
                 outlier.size = outlier.size),
    facet_rep_wrap(~PFT, scales = "free", ncol = 2),
    scale_fill_manual(values = cols_graze, 
                      name = "Grazing Treatment"),
    # so just display the RCP
    scale_x_discrete(labels = years2lab),
    theme(legend.position = c(0.75, 0.15)),
    geom_vline(xintercept = 2.5, linetype = 2),
    labs(y = lab_bio2,
         x = lab_yrs)
    
  )
}

jpeg("figures/biomass/pub_qual/pft5_bio_diff_boxplot_c4on.jpeg",
     res = 600, height = 8, width = 5, units = "in")

pft5_bio_d2 %>% 
  filter(c4 == "c4on") %>% 
  ggplot(aes(id2, bio_diff, fill = graze)) +
  box2()

dev.off()

jpeg("figures/biomass/pub_qual/pft5_bio_diff_boxplot_c4off.jpeg",
     res = 600, height = 8, width = 5, units = "in")

pft5_bio_d2 %>% 
  filter(c4 == "c4off") %>% 
  ggplot(aes(id2, bio_diff, fill = graze)) +
  box2()

dev.off()

# * change relative to reference graze ------------------------------------
# boxplots showing change in biomass relative to current time period and
# given grazing intensity


# ** boxplot ---------------------------------------------------------------

pdf("figures/biomass/pft5_bio_diff_gref_boxplots.pdf", 
    height = 6.5, width = 5)

# 'loop' over referenc class
map2(pft5_d_grefs, names(pft5_d_grefs), function(df, ref_graze){
  
  # 'loop' over c4 levs
  map(levs_c4, function(lev_c4){
  
    df <- df %>% 
    # just c4 on for now, adjust the line and change file name to run for c4off
    filter(c4 == lev_c4)

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
    scale_fill_graze() +
    # so just display the RCP
    scale_x_discrete(labels = years2lab) +
    theme(legend.position = c(0.75, 0.15)) +
    geom_vline(xintercept = c(1.5, 3.5), linetype = 2) +
    labs(y = lab_bio2,
         x = lab_yrs,
         subtitle = paste("Change in biomass relative to", tolower(ref_graze), 
                          "grazing \n under current conditions"),
         caption = c4on_off_lab(lev_c4))
g
})
})
dev.off()


# ** scatterplot (vs climate) ----------------------------------------------

pdf("figures/biomass/bio-diff_vs_climate_v1.pdf",
    width = 6, height = 5)

# bio-diff vs MAP and MAT, for for ref class of light grazing, for
# each RCP/time period and PFT

pmap(levs_pft_c4, function(pft, levs_c4) {
  # using Light grazing as reference class
  g <- pft5_d_grefs[["Light"]] %>% 
    filter(PFT == pft, c4 == levs_c4) %>% 
    ggplot(aes(y = bio_diff, color = graze)) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_rep_wrap(~RCP + years) +
    scale_color_graze() +
    labs(y = lab_bio2,
         caption = paste(c4on_off_lab(levs_c4),
           "\nReference class is light grazing under current conditions"),
         subtitle = paste("Change in", pft, "biomass")) +
    theme(legend.position = c(0.85, 0.15),
          axis.text = element_text(size = 7))
  
  climate_scatter(g)
  
})

dev.off()


# fire --------------------------------------------------------------------

# boxplot of fire return interval by RCP, and grazing intensity

map(levs_c4, function(lev_c4) {
  g1 <- fire1 %>% 
    filter(c4 == lev_c4) %>% 
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
    filter(c4 == lev_c4) %>% 
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
    scale_fill_graze() +
    # so just display the RCP
    scale_x_discrete(labels = years2lab) +
    theme(legend.position = "top") +
    geom_vline(xintercept = 2.5, linetype = 2) +
    labs(x = lab_yrs,
         y = lab_fire1,
         caption = "Outliers not shown in either panel")
  
  jpeg(paste0("figures/fire/fire_return_boxplots_", lev_c4, ".jpeg"),
       height = 8, width = 5, res = 600, units = "in")
  gridExtra::grid.arrange(
    g1, g2
  )
  dev.off()
})
