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
line_loc2 <- c(1.5, 3.5) # for figures w/ 2 vertical lines
outlier.size = 0.5

# colors
scale_color_graze <- function() {
  scale_color_manual(values = cols_graze, name = "Grazing")
}

scale_fill_graze <- function() {
  scale_fill_manual(values = cols_graze, name = "Grazing")
}

# width of figs
wfig_box1 <- 9 # width of boxplots in inches

# location of legends in boxplots
legend_pos_box1 <- "top"

# number of columns of panels in boxplots
ncol_box <- 3

# * vectors/dfs for 'looping' ---------------------------------------------
# for looping
levs_pft <- levels(pft5_bio2$PFT)

levs_pft_c4 <- expand_grid(pft = factor(levs_pft, levs_pft), 
                           levs_c4 = unique(pft5_bio2$c4)) %>% 
  arrange(pft, desc(levs_c4))

levs_c4 <- unique(levs_pft_c4$levs_c4) %>% 
  sort(., decreasing = TRUE)

levs_grefs_c4 <- expand_grid(ref_graze = names(pft5_d_grefs), 
                             levs_c4 = levs_c4)

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

# for adding a second y axis to effect size figures, that shows % change
add_sec_axis <- function() {
  scale_y_continuous(sec.axis = sec_axis(trans = es2pchange, 
                                         name = "% Change")) 
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
    facet_rep_wrap(~ PFT, scales = "free", ncol = ncol_box), 
    scale_fill_manual(values = cols_rcp, name = "Scenario"),
    scale_x_discrete(labels = years2lab),
    geom_vline(xintercept = line_loc, linetype = 2),
    theme(legend.position = legend_pos_box1,
          axis.text = element_text(size = 7)),
    labs(x = lab_yrs,
         y = lab_bio0)
  )
  }


jpeg("figures/biomass/pub_qual/pft5_bio_boxplot_c4on.jpeg",
     res = 600, height = 8, width = wfig_box1, units = "in")

pft5_bio2 %>% 
  filter(c4 == "c4on") %>% 
  # so boxplot doesn't show sites with 0 biomass
  mutate(biomass = ifelse(biomass == 0, NA, biomass)) %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  box1()

dev.off()

jpeg("figures/biomass/pub_qual/pft5_bio_boxplot_c4off.jpeg",
     res = 600, height = 8, width = wfig_box1, units = "in")

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

# boxplot of change in biomass (scaled percent and effect size), 
# for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period

box2 <- function(axis_data, var = "bio_diff", mult = 0.05,
                 subtitle = "Comparing to current conditions within a grazing level",
                 xintercept = 2.5,
                 repeat.tick.labels = FALSE){
  list(
    geom_boxplot(position = "dodge",
                 outlier.size = outlier.size),
    facet_rep_wrap(~PFT, scales = "free_y", ncol = ncol_box,
                   repeat.tick.labels = repeat.tick.labels),
    scale_fill_graze(),
    # so just display the RCP
    scale_x_discrete(labels = years2lab),
    theme(legend.position = legend_pos_box1),
    geom_vline(xintercept = xintercept, linetype = 2),
    labs(x = lab_yrs, subtitle = subtitle),
    # text and empty points based on a different dataframe, so that
    # axis limits amongst multiple figures can be the same
    geom_text(data = box_anno(axis_data, var = var, 
                               group_by = c("PFT", "RCP"),
                               id = "id2", mult = mult),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5),
    geom_point(data = axis_data, color = NA)
    
  )
}


pdf("figures/biomass/pft5_bio_diff_boxplots_v1.pdf",
    height = 8, width = wfig_box1)

# % change
map(levs_c4, function(lev_c4) {

  pft5_bio_d2 %>% 
    filter(c4 == lev_c4) %>% 
    ggplot(aes(id2, bio_diff, fill = graze)) +
    box2(axis_data = pft5_bio_d2) +
    labs(y = lab_bio2,
         caption = c4on_off_lab(lev_c4))

})

# effect size
map(levs_c4, function(lev_c4) {

  pft5_bio_es1 %>% 
    filter(c4 == lev_c4) %>% 
    ggplot(aes(id2, bio_es, fill = graze)) +
    box2(axis_data = pft5_bio_es1, 
         var = "bio_es", repeat.tick.labels = "y") +
    add_sec_axis() +
    labs(y = lab_es0,
         caption = c4on_off_lab(lev_c4))
  
})

dev.off()

# * change relative to reference graze ------------------------------------


# ** boxplot ---------------------------------------------------------------
# boxplots showing change in biomass (and effect size) relative to current time 
# period and a given grazing intensity ('reference class')

box3 <- function(axis_data, var = "bio_diff") {
  list(
      geom_hline(yintercept = 0, alpha = 0.3, linetype = 1),
      geom_boxplot(position = position_dodge(preserve = "single"),
                   outlier.size = outlier.size),
      facet_rep_wrap(~PFT, scales = "free", ncol = ncol_box),
      scale_fill_graze(),
      # so just display the RCP
      scale_x_discrete(labels = years2lab),
      theme(legend.position = legend_pos_box1),
      geom_vline(xintercept = line_loc2, linetype = 2),
      labs(x = lab_yrs),
      geom_text(data = ~box_anno(axis_data, var = var, 
                                 group_by = c("PFT", "RCP"),
                                 id = "id2"),
                aes(x, y, label = RCP, fill = NULL),
                size = 2.5),
      geom_point(data = axis_data, color = NA)
  )
}

pdf("figures/biomass/pft5_bio_diff_gref_boxplots_v1.pdf", 
    height = 6.5, width = wfig_box1)

# % change
# 'loop' over reference class (reference grazing level) and c4 on or off
pmap(levs_grefs_c4, function(ref_graze, levs_c4){
  
  df0 <- pft5_d_grefs[[ref_graze]]
  
  df <- df0 %>% 
    filter(c4 == levs_c4) 

  # at this point it is not possible to show the unused grazing level in ggplot
  # for the current rcp (https://github.com/tidyverse/ggplot2/issues/3345)
  ggplot(df, aes(id2, bio_diff, fill = graze)) +
    # keeping axes the same between the c4 on and off figures
    box3(axis_data = df0) +
    labs(y = lab_bio2,
         subtitle = paste("Change in biomass relative to", tolower(ref_graze), 
                          "grazing \n under current conditions"),
         caption = c4on_off_lab(levs_c4))
})

# effect size
pmap(levs_grefs_c4, function(ref_graze, levs_c4){
  
  df0 <- pft5_es_grefs[[ref_graze]]
  
  df <- df0 %>% 
    filter(c4 == levs_c4) 
  
  ggplot(df, aes(id2, bio_es, fill = graze)) +
    box3(axis_data = df0, var = "bio_es") +
    add_sec_axis() +
    labs(y = lab_es0,
         subtitle = paste("Change in biomass relative to", tolower(ref_graze), 
                          "grazing \n under current conditions"),
         caption = c4on_off_lab(levs_c4))
})


dev.off()


# ** scatterplot (vs climate) ----------------------------------------------

scatter1 <- function(pft, levs_c4, axis_data) {
  list(
    geom_hline(yintercept = 0, linetype = 2),
    facet_rep_wrap(~RCP + years),
    scale_color_graze(),
    labs(caption = paste(c4on_off_lab(levs_c4),
                         "\nReference class is light grazing under current",
                         "conditions"),
         subtitle = paste("Change in", pft, "biomass")),
    theme(legend.position = c(0.85, 0.15),
          axis.text = element_text(size = 7)),
    # for setting axis limits based on bigger data frame, so multiple
    # figs can have same y axis lims
    geom_point(data = axis_data,
               x = NA)
  )
}

pdf("figures/biomass/bio-diff_vs_climate_v1.pdf",
    width = 6, height = 5)

# bio-diff and effect size vs MAP and MAT, for for ref class of light grazing, for
# each RCP/time period and PFT

# %change
l1 <- pmap(levs_pft_c4, function(pft, levs_c4) {
  # using Light grazing as reference class
  df <- pft5_d_grefs[["Light"]] %>% 
    filter(PFT == pft)
  
  g <- df %>% 
    filter(c4 == levs_c4) %>% 
    ggplot(aes(y = bio_diff, color = graze)) +
    scatter1(pft, levs_c4, axis_data = df) +
    labs(y = lab_bio2)
  
  climate_scatter(g)
  
})


# effect size
l2 <- pmap(levs_pft_c4, function(pft, levs_c4) {
  # using Light grazing as reference class
  df <- pft5_es_grefs[["Light"]] %>% 
    filter(PFT == pft)
  
  g <- df %>% 
    filter(c4 == levs_c4) %>% 
    ggplot(aes(y = bio_es, color = graze)) +
    scatter1(pft, levs_c4, axis_data = df) +
    labs(y = lab_es0)
  
  climate_scatter(g)
  
})

map_depth(c(l1, l2), .depth = 1, .f = `[`, "MAT") # all MAT figs
map_depth(c(l1, l2), .depth = 1, .f = `[`, "MAP") # all MAP figs

dev.off()


# * change w/ in a scenario ----------------------------------------------

# change in biomass going from light grazing, to some other grazing level,
# for a given scenario (e.g. light grazing RCP 8.5 end century to heavy graze
# RCP 8.5 end century)


# ** boxplots -------------------------------------------------------------

pdf("figures/biomass/pft5_bio_diff_wgcm_boxplots_v1.pdf", 
    height = 6.5, width = wfig_box1)

# % change
map(levs_c4, function(x){
  
  pft5_d_wgcm %>% 
    filter(c4 == x) %>% 
    ggplot(aes(id2, bio_diff, fill = graze)) +
    box2(axis_data = pft5_d_wgcm, xintercept = line_loc2,
         subtitle = "Reference group is light grazing for the given climate scenario") +
    labs(y = lab_bio2,
         caption = c4on_off_lab(x))
  
})

map(levs_c4, function(x){
 
   pft5_es_wgcm %>% 
    filter(c4 == x) %>% 
    ggplot(aes(id2, bio_es, fill = graze)) +
    box2(axis_data = pft5_es_wgcm, var = "bio_es", xintercept = line_loc2,
         subtitle = "Reference group is light grazing for the given climate scenario") +
    add_sec_axis() +
    labs(y = lab_es0,
         caption = c4on_off_lab(x))
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
