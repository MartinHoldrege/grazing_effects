# Martin Holdrege

# Script started Jan 31, 2022

# Purpose: 
# Create exploratory figures of the number of individual plants by 
# by plant functional type ('indivs' column)
# to assess how STEPWAT2 simulations with c4 expansion on vs off
# differ in number of individuals
# One main question is to understand whether biomass is being constrained
# by the number of individuals, and whether there is much
# variability in the number of individuals between sites etc.
# Second objective is to see how the number of individuals is affected
# by climate scenarios and grazing treatments

# Figure ideas (delete once created):

# dependencies ------------------------------------------------------------

library(lemon) # for face_rep_wrap()
library(tidyverse)
path_rds <- 'data_processed/site_means/summarize_bio.RDS'

if(!file.exists(path_rds)) {
  source("scripts/02_summarize_bio.R") # where needed dataframes are created
}

source("src/fig_params.R")
source("src/fig_functions.R")
source('src/general_functions.R')

# params ------------------------------------------------------------------
# misc. parameters and settings for figure creation

theme_set(theme_custom1())
theme_update(strip.background = element_blank())
line_loc2 <- c(1.5, 3.5) # for figures w/ 2 vertical lines

# width of figs
wfig_box1 <- 9 # width of boxplots in inches

# number of columns of panels in boxplots
ncol_box <- 3


# dataframes etc --------------------------------------------------------------

l <- readRDS(path_rds)

pft5_bio2 <- l$pft5_bio2
pft5_bio_d2 <- l$pft5_bio_d2
pft5_bio_es1 <- l$pft5_bio_es1
pft5_d_wgcm <- l$pft5_d_wgcm
pft5_es_wgcm <- l$pft5_es_wgcm
runs_graze <- l$runs_graze
pft5_d_grefs <- l$pft5_d_grefs
group_cols <- l$group_cols

# * vectors/dfs for 'looping' ---------------------------------------------
# for looping
levs_pft0 <- levels(pft5_bio2$PFT) # all pft levels in the main datafile

# only include levels actually in the dataset
levs_pft <- levs_pft0[levs_pft0 %in% unique(pft5_bio2$PFT)]

levs_pft_run <- expand_grid(pft = factor(levs_pft, levs_pft0), 
                            run = runs_graze) %>% 
  arrange(pft, run)


levs_grefs_run <- expand_grid(ref_graze = names(pft5_d_grefs), 
                              run = runs_graze)


# list of pfts
levs_pft_l <- list(pft5 = pft5_factor(x = NULL, return_levels = TRUE)) # the main 5

# other pfts
levs_pft_l$other = levs_pft[!levs_pft %in% levs_pft_l$pft5]

# absolute indivs ---------------------------------------------------------

# figures of individuals (absolute numbers not % change etc), by grazing and
# climate scenario

# * boxplot ---------------------------------------------------------------

# c4 on and off in side by side panels 
# for each pft
pdf("figures/indivs/indivs_boxplot_by_run_v2.pdf",
    height = 9, width = 6)

map(levs_pft_l, function(x) {
  df <- pft5_bio2 %>% 
    filter(PFT %in% x) %>% 
    mutate(indivs = ifelse(indivs == 0, NA, indivs))

  ggplot(df, aes(x = id, y = indivs, fill = RCP)) +
    box1(var = "indivs", 
         add_facet_wrap = FALSE,
         y = lab_indivs0) +
    facet_rep_grid(PFT ~ run, scales = "free_y") +
    labs(caption = "Note: Data points w/ 0 individuals not shown")

})

dev.off()

# * scatterplot (vs climate) ----------------------------------------------

# each combination PFt and C4 on or off is plotted on seperate pages
# facets are RCP/year combinations, colors are grazing
pdf("figures/indivs/indivs_vs_climate_v3.pdf",
    width = 6, height = 5)

# here pmap is working by taking the column names (i.e. iterating over rows)
# individuals vs map or MAT by grazing level, facets are climate
# scenario
pmap(levs_pft_run, function(pft, run) {
  g <- pft5_bio2 %>% 
    filter(PFT == pft, run == run) %>% 
    ggplot(aes(y = indivs, color = graze)) +

    facet_rep_wrap(~RCP + years) +
    scale_color_graze() +
    labs(y = lab_indivs0,
         subtitle = paste(pft, "individuals"),
         caption = run) +
    theme(legend.position = c(0.85, 0.15),
          axis.text = element_text(size = 7)) 
  
  climate_scatter(g)
})

# facets are grazing level, colors are RCP
pmap(levs_pft_run, function(pft, run) {
  g <- pft5_bio2 %>% 
    filter(PFT == pft, run == run,
           years != "2030-2060") %>% 
    ggplot(aes(y = indivs, color = RCP)) +
    facet_rep_wrap(~graze) +
    scale_color_manual(values = cols_rcp) +
    labs(y = lab_indivs0,
         subtitle = paste(pft, "individuals"),
         caption = paste(run,
                         "\nOnly showing 2070-2100 for RCP 4.5 & 8.5")) +
    theme(legend.position = "top",
          axis.text = element_text(size = 7))
  
  climate_scatter(g)
  
})

dev.off()

# c4 on vs off indivs ----------------------------------------------------

# code for figures showing c4on vs off (ie scatterplots comparing runs)
# are available in commits pre 2025. 

# indivs change -------------------------------------------------------
# response variable is change in individuals

# * change relative to same graze ------------------------------------

# boxplot of change in individuals (scaled percent and effect size), 
# for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period, relative to current
# conditions for that same grazing level (i.e. this is the climate effect)


pdf("figures/indivs/pft5_indivs_diff_boxplots_v4.pdf",
    height = 8, width = wfig_box1)

# % change
map(runs_graze, function(run) {
  
  pft5_bio_d2 %>% 
    filter(run == !!run) %>% 
    ggplot(aes(id2, indivs_diff, fill = graze)) +
    box2(axis_data = pft5_bio_d2,
         var = "indivs_diff") +
    labs(y = lab_indivs2,
         caption = run)
  
})

# effect size (outliers not shown)
map(runs_graze, function(run) {
  
  # compute boxplot statistics, so can plot without outliers (many extreme
  # outliers)
  df <- pft5_bio_es1 %>% 
    # group_cols vector define in 02_summarize_bio.R
    group_by(across(all_of(c(group_cols[!group_cols %in% c("site", "GCM")],
                             "id2")))) %>% 
    compute_boxplot_stats(var = "indivs_es")
  
  df %>% 
    filter(run == !!run) %>% 
    ggplot(aes(x = id2, fill = graze)) +
    box2(axis_data =  boxplot_stats_long(df), # so runs have same limits
         var = "y", repeat.tick.labels = "y",
         box_identity = TRUE) +
    add_sec_axis() +
    labs(y = lab_es0_indivs,
         caption = paste("Outliers not shown\n", run))
  
})

dev.off()

# * wgcm ----------------------------------------------

# change in num of indivs going from light grazing, to some other grazing level,
# for a given scenario, ie. change within a gcm/climate scenario
# (e.g. light grazing RCP 8.5 end century to heavy graze, RCP 8.5 end century)


# ** boxplots -------------------------------------------------------------

pdf("figures/indivs/pft5_indivs-diff_wgcm_boxplots_v5.pdf", 
    height = 8, width = wfig_box1)

# % change
map(runs_graze, function(run){
  
  pft5_d_wgcm %>% 
    filter(run == !!run) %>% 
    ggplot(aes(id2, indivs_diff, fill = graze)) +
    box2(axis_data = pft5_d_wgcm, 
         var = "indivs_diff",
         xintercept = line_loc2,
         subtitle = "Reference group is light grazing for the given climate scenario") +
    labs(y = lab_indivs2,
         caption = run)
  
})

# effect size
map(runs_graze, function(run){
  
  out <- list()
  
  # w/ outliers
  out[[1]] <- pft5_es_wgcm %>% 
    filter(run == !!run) %>% 
    ggplot(aes(id2, indivs_es, fill = graze)) +
    box2(axis_data = pft5_es_wgcm, 
         var = "indivs_es", xintercept = line_loc2,
         subtitle = "Reference group is light grazing for the given climate scenario",
         repeat.tick.labels = TRUE) +
    add_sec_axis() +
    labs(y = lab_es0_indivs,
         caption = paste0("All data shown\n", run)) +
    theme(axis.text.y.right = element_text(size = 5)) 
  
  # boxplot outliers removed
  df <-  pft5_es_wgcm %>% 
    group_by(id2, PFT, RCP, graze, years, run) %>% 
    compute_boxplot_stats(var = "indivs_es")
  
  g <- df %>% 
    filter(run == !!run) %>% 
    ggplot(aes(x = id2, fill = graze)) +
    add_sec_axis() +
    labs(y = lab_es0_indivs,
         caption = paste0("Outliers removed\n", run))
  
  out[[2]] <- g +
    box2(axis_data = boxplot_stats_long(df), var = "y", xintercept = line_loc2,
         box_identity = TRUE, 
         subtitle = "Reference group is light grazing for the given climate scenario",
         repeat.tick.labels = TRUE
    )  
  
  # same figure but with fixed scales
  out[[3]] <- g +
    box2(axis_data = boxplot_stats_long(df), var = "y", xintercept = line_loc2,
         box_identity = TRUE, 
         subtitle = 
           paste("Y axis Scales fixed.", 
                 "Reference group is light grazing for the given climate scenario"),
         repeat.tick.labels = TRUE,
         scales = "fixed",
         anno_same_across_panels = TRUE)
  
  out
})

dev.off()

