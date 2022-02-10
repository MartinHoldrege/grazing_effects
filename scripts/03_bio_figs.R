# Martin Holdrege

# Script started 1/6/2022

# Purpose--to create figures of biomass data from simulations
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
line_loc2 <- c(1.5, 3.5) # for figures w/ 2 vertical lines

# width of figs
wfig_box1 <- 9 # width of boxplots in inches
wfig_box2 <- 6 # width of boxplots that show just 5 pft groups
hfig_box2 <- 5 # height of 5 panel boxplots

# number of columns of panels in boxplots
ncol_box <- 3

# * vectors/dfs for 'looping' ---------------------------------------------
# for looping
levs_pft <- levels(pft5_bio2$PFT) # all pft levels in the main datafile

# for when need to 'loop' over both pft4 and c4 on and off
levs_pft_c4 <- expand_grid(pft = factor(levs_pft, levs_pft), 
                           levs_c4 = unique(pft5_bio2$c4)) %>% 
  arrange(pft, desc(levs_c4))

levs_c4 <- unique(levs_pft_c4$levs_c4) %>% 
  sort(., decreasing = TRUE)

levs_grefs_c4 <- expand_grid(ref_graze = names(pft5_d_grefs), 
                             levs_c4 = levs_c4)

# list of pfts
levs_pft_l <- list(pft5 = pft5_factor(x = NULL, return_levels = TRUE)) # the main 5

# other pfts
levs_pft_l$other = levs_pft[!levs_pft %in% levs_pft_l$pft5]


# functions ---------------------------------------------------------------
# If function needs to be used another script move to fig_functions.R

# for boxplots showing only the main 5 PFTs
theme_box_pft5 <- function() {
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.6),
        legend.position = c(0.85, 0.2))
 }

# prep data ---------------------------------------------------------------

# data frame to be used for boxplots
pft5_bio_b <- pft5_bio2 %>% 
  # consider removing 0s so boxplot doesn't show sites with 0 biomass
 mutate(biomass = ifelse(biomass == 0, NA, biomass))

# absolute biomass ------------------------------------------------------
# This figure meant to be analogous to M.E.'s thesis figure 9.
# biomass by pft, rcp, time period, and grazing intensity


# * boxplot ---------------------------------------------------------------

jpeg("figures/biomass/pub_qual/bio-boxplot_all_c4on.jpeg",
     res = 600, height = 8, width = wfig_box1, units = "in")

pft5_bio_b %>% 
  filter(c4 == "c4on") %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  box1()

dev.off()

jpeg("figures/biomass/pub_qual/bio-boxplot_all_c4off.jpeg",
     res = 600, height = 8, width = wfig_box1, units = "in")

pft5_bio_b %>% 
  filter(c4 == "c4off") %>% 
  ggplot(aes(x = id, y = biomass, fill = RCP)) +
  box1()

dev.off()

# c4 on and off in side by side panels 
# for each pft
pdf("figures/biomass/bio_boxplot_c4off_and_on.pdf",
     height = 9, width = 6)

map(levs_pft_l, function(x) {
  df <- pft5_bio2 %>% 
    filter(PFT %in% x)
  out <- list()
  out[[1]] <- df %>%  # removing 0s
    mutate(biomass = ifelse(biomass == 0, NA, biomass)) %>% 
    ggplot(aes(x = id, y = biomass, fill = RCP)) +
    box1(add_facet_wrap = FALSE) +
    facet_rep_grid(PFT ~ c4, scales = "free_y") +
    labs(caption = "Note: 0 biomass values not shown")
  
  # showing 0s
  out[[2]] <- ggplot(df, aes(x = id, y = biomass, fill = RCP)) +
    box1(add_facet_wrap = FALSE) +
    facet_rep_grid(PFT ~ c4, scales = "free_y") +
    labs(caption = "Note: 0 biomass values are shown")
  out
})

dev.off()

# ** fewer scenarios ------------------------------------------------------

# showing the biomass boxplot but only for the main PFTs, and RCP mid-century
# for possible use in manuscript

df <- pft5_bio_b %>% 
  filter_rcp_c4(PFT = TRUE)

n_zero <- pft5_bio2 %>% 
  filter_rcp_c4(PFT = TRUE) %>% 
  group_by(RCP, PFT, graze) %>% 
  summarize(n_zero = sum(biomass == 0, na.rm =TRUE),
            .groups = 'drop_last')

# check
if (any(summarize(n_zero, lu = lu(n_zero))$lu != 1)) {
  warning("different number of sites with 0 biomass between grazing trmts")
} 

# text to add to figure regarding how many sites had 0 biomass
n_zero2 <- n_zero %>% 
  summarise(n_zero = unique(n_zero),
            .groups = 'drop') %>% 
  mutate(string = paste0(n_zero, " sites (", n_zero/200*100, "%)\n had 0 biomass"),
         string = ifelse(n_zero == 0, "", string),
         # singular
         string = ifelse(n_zero == 1, str_replace(string, "sites", "site"),
                         string),
         graze = 'Light',
         y = ifelse(PFT == "Cheatgrass", -6, 0) #  where to put label
         )
n_zero2

g <- ggplot(df, aes(x = RCP, y = biomass, fill = graze)) +
  facet_rep_wrap(~ PFT, scales = "free", ncol = 3) +
  # currently seems like I need to restrict the named color vector to the names
  # being plotted, see https://github.com/tidyverse/ggplot2/pull/4619
  scale_fill_graze() +
  labs(x = lab_rcp,
       y = lab_bio0) +
  theme_box_pft5() +
#  ylim(-5, NA) +
  geom_text(data = n_zero2, aes(y = y, label = string),
            size = 2,vjust = 'inward',
            lineheight = 0.7)


# boxplot
jpeg("figures/biomass/pub_qual/bio-boxplot_pft5_rcp8.5_c4on.jpeg",
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")
g +  geom_boxplot(outlier.size = outlier.size) +
  # for some reason geom_vline needs to come after boxplot, or error thrown
  geom_vline(xintercept = 1.5, linetype = 2)


dev.off()

jpeg("figures/biomass/pub_qual/bio-violin_pft5_rcp8.5_c4on.jpeg",
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")

g + geom_violin() +
  geom_vline(xintercept = 1.5, linetype = 2)

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


# * cheat light vs heavy graze --------------------------------------------
# comparing cheatgrass light vs heavy grazing, to help determine why some
# sites have such large % changes in cheatgrass (it appears it is the
# group of sites with low cheat biomass)

cheat_df <- pft5_bio2 %>% 
  filter(PFT == "Cheatgrass", RCP %in% c("Current", "RCP8.5")) %>% 
  select(c4, years, RCP, graze, site, biomass)

pdf("figures/biomass/light_vs_heavy_graze_scatterplot.pdf",
    width = 6, height = 6)
cheat_df %>% 
  filter(graze == "Light") %>% 
  left_join(filter(cheat_df, graze == "Heavy"),
            by = c("c4", "RCP", "years", "site"),
            suffix = c("_light", "_heavy")) %>% 
  ggplot(aes(biomass_light, biomass_heavy)) +
  geom_point() +
  facet_rep_wrap(~RCP + years, ncol = 2) +
  geom_abline(slope = 1) +
  labs(x = "Biomass (light grazing)",
       y = "Biomass (heavy grazing)",
       subtitle = "Cheatgrass biomass under light vs heavy grazing")
dev.off()
# c4 on vs off biomass ----------------------------------------------------

pdf("figures/biomass/bio_c4on_vs_off_v1.pdf",
    width = 9, height = 7)


# * absolute biomass ------------------------------------------------------

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


# * on vs off %change ----------------------------------------------------
# boxplot of % change of going from c4 on to c4 off
c4on_v_off_diff %>% 
  #filter(PFT != "C4Pgrass") %>% 
  group_by(years, RCP, graze, PFT, id) %>% 
  compute_boxplot_stats(var = "bio_es") %>% 
  ggplot(aes(x = id, fill = RCP)) +
  geom_text(data = ~box_anno(boxplot_stats_long(.), 
                             var = "y", group_by = c("PFT", "graze"),
                             mult = 0.05),
            aes(x, y, label = graze, fill = NULL), 
            size = 2.5) +
  geom_boxplot_identity() +
  scale_fill_manual(values = cols_rcp, name = "Scenario") +
  scale_x_discrete(labels = id2year) +
  geom_vline(xintercept = line_loc, linetype = 2) +
  theme(legend.position = legend_pos_box1,
        axis.text = element_text(size = 7)) +
  facet_rep_wrap(~ PFT, scales = "free", ncol = ncol_box,
                 repeat.tick.labels = TRUE) +
  add_sec_axis() +
  labs(x = lab_yrs,
       y = lab_es_on_off,
       caption = "Outliers not shown. 
       Not showing points where biomass was 0 for on or off",
       subtitle = "Change in biomass when C4 expansion is turned off 
       in a given climate scenario and grazing level") 


dev.off()

# C3Pgrass/Pgrass ---------------------------------------------------------

g <- ggplot(C3_Pgrass_ratio, aes(id, C3_Pgrass_ratio, fill = RCP)) +
  box1(add_facet_wrap = FALSE, var = "C3_Pgrass_ratio",
       group_by = c("c4", "graze"),
       y = lab_C3_Pgrass_ratio) +
  facet_rep_wrap(~c4) +
  labs(subtitle = "Ratio of C3Pgrass biomass to total Pgrass biomass",
       caption = "Separately showing data from simulations with C4 expansion on and off") 

ggsave("figures/biomass/C3_Pgrass_ratio.jpeg", g, 
       width = 6, height = 4)

# biomass change -------------------------------------------------------
# response variable is change in biomass

# * change relative to same graze ------------------------------------

# boxplot of change in biomass (scaled percent and effect size), 
# for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period, relative to current
# conditions for that same grazing level (i.e. this is the climate effect)


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
    scatter_light(pft, levs_c4, axis_data = df) +
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
    scatter_light(pft, levs_c4, axis_data = df) +
    labs(y = lab_es0)
  
  climate_scatter(g)
  
})

map_depth(c(l1, l2), .depth = 1, .f = `[`, "MAT") # all MAT figs
map_depth(c(l1, l2), .depth = 1, .f = `[`, "MAP") # all MAP figs

dev.off()


# * wgcm ----------------------------------------------

# change in biomass going from light grazing, to some other grazing level,
# for a given scenario, ie. change within a gcm/climate scenario
# (e.g. light grazing RCP 8.5 end century to heavy graze, RCP 8.5 end century)


# ** boxplots -------------------------------------------------------------

pdf("figures/biomass/pft5_bio_diff_wgcm_boxplots_v1.pdf", 
    height = 8, width = wfig_box1)

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

# effect size
map(levs_c4, function(x){
 
  out <- list()
  # w/outliers
   out[[1]] <- pft5_es_wgcm %>% 
    filter(c4 == x) %>% 
    ggplot(aes(id2, bio_es, fill = graze)) +
    box2(axis_data = pft5_es_wgcm, 
         var = "bio_es", xintercept = line_loc2,
         subtitle = "Reference group is light grazing for the given climate scenario",
         repeat.tick.labels = TRUE) +
    add_sec_axis() +
    labs(y = lab_es0,
         caption = paste0("All data shown\n", c4on_off_lab(x))) +
    theme(axis.text.y.right = element_text(size = 5)) 
  
  # boxplot outliers removed
   df <-  pft5_es_wgcm %>% 
     group_by(id2, PFT, RCP, graze, years, c4) %>% 
     compute_boxplot_stats(var = "bio_es")
  
   g <- df %>% 
     filter(c4 == x) %>% 
     ggplot(aes(x = id2, fill = graze)) +
     add_sec_axis() +
     labs(y = lab_es0,
          caption = paste0("Outliers removed\n", c4on_off_lab(x)))
   
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
         subtitle = "Reference group is light grazing for the given climate scenario",
         repeat.tick.labels = TRUE,
         scales = "fixed",
         anno_same_across_panels = TRUE)
     
  out
})

dev.off()


# ***fewer scenarios ------------------------------------------------------

# effect size boxplot/violin plot for main PFTs, and only RCP8.5 mid-century,
# (i.e. closer to publication quality)


ylim <- c(-1.9, 1) # limits for cheatgrass


df <- pft5_es_wgcm %>% 
  filter_rcp_c4(PFT = TRUE)

outliers <- df%>% 
  filter(bio_es < ylim[1])

if(!all(outliers$PFT == "Cheatgrass")) {
  stop("Other PFTs also have outliers")
} else {
  outliers <- outliers %>% 
    group_by(graze, RCP, PFT) %>% 
    summarize(n = n(),
              .groups = 'drop') %>% 
    mutate(n = paste0(n, "*"))
}

# base of the plot
g <- ggplot(df, aes(x = RCP, y = bio_es, fill = graze)) +
  geom_blank() + # added so that geom_vline doesn't throw an error
  geom_vline(xintercept = 1.5, linetype = 2) +
  facet_rep_wrap(~ PFT, scales = "free", ncol = 3,
                 repeat.tick.labels = TRUE) +
  # setting axis limits separately for the cheatgrass panel
  ggh4x::facetted_pos_scales(y = list(
    PFT == "Cheatgrass" ~ scale_y_continuous(
      limits = ylim, #
      sec.axis = sec_axis(trans = es2pchange,name = "% Change",
                          breaks = c(-75, 0, 75, 150)))
  ))+
  scale_fill_graze(include_light = FALSE) +
  labs(x = lab_rcp,
       y = lab_es1) +
  theme_box_pft5()+
  add_sec_axis() +
  # so non cheatgrass panels all have the same limits 
  expand_limits(y = c(-1.1, 0.35)) +
  # printing how many outliers not shown. 
  geom_text(data = outliers, aes(x = as.numeric(RCP) + as.numeric(graze)/3.7 -0.5,
                                 y = -2, label = n, color = graze),
            size = 3, vjust = 'inward') +
  scale_color_manual(values = cols_graze) +
  guides(color = 'none')

# boxplot
jpeg("figures/biomass/pub_qual/bio-diff_boxplot_pft5_rcp8.5_c4on.jpeg",
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")
g +
  geom_boxplot(outlier.size = outlier.size) 

dev.off()

# violin plot
jpeg("figures/biomass/pub_qual/bio-diff_violin_pft5_rcp8.5_c4on.jpeg",
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")

# applying filter otherwise violin goes off the edge of the page
g +
  geom_violin(data = df %>% filter(bio_es > ylim[1])) 

dev.off()

# % below threshold ------------------------------------------------------

ref_threshold2 <- ref_threshold %>% 
  mutate(string = paste(round(threshold, 0), "g/m^2"))

pdf("figures/biomass/threshold_dotplots_v1.pdf", 
    height = 6.5, width = wfig_box1)

map(levs_c4, function(lev_c4) {
  
threshold1 %>% 
  filter(c4 == lev_c4) %>% 
  ggplot(aes(x = id, y = pcent_above)) +
  geom_text(data = ~box_anno(., var = "pcent_above",
                             group_by = c("PFT", "graze"),
                             y = 110),
            aes(x, y, label = graze, fill = NULL),
            size = 2.5) +
  # printing the value of the threshold on each panel
  geom_text(data = filter(ref_threshold2, c4 == lev_c4),
            aes(x = 3, y = 10, label = string),
                size = 2.5) +
  # line at the threshold
  geom_hline(yintercept = (1 - pcent)*100, alpha = 0.5) +
  geom_point(aes(color = RCP)) +
  facet_rep_wrap(~ PFT, ncol = ncol_box) +
  scale_color_manual(values = cols_rcp, name = "Scenario") +
  scale_x_discrete(labels = id2year) +
  geom_vline(xintercept = line_loc, linetype = 2) +
  theme(legend.position = legend_pos_box1,
        axis.text = element_text(size = 7)) +
  labs(x = lab_yrs,
       y = lab_below0,
       subtitle = "% of sites with biomass above current light grazing 5th percentile",
       caption = paste0(c4on_off_lab(lev_c4),
                       "\nValue of threshold provided in each panel"))
  
})

dev.off()


# composition -------------------------------------------------------------
# stacked barcharts showing the composition of the 5 main PFTs.
# this is based on the mean biomass across sites for each of the PFTs

g <- ggplot(comp2, aes(x = id, fill = PFT)) +
  scale_x_discrete(labels = id2graze) +
  labs(x = lab_graze) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) +
  geom_vline(xintercept = 4.5, linetype = 2) 

# absolute biomass (stacked)
g1 <- g +
  geom_bar(aes(y = biomass), position = 'stack', stat = 'identity') +
  geom_text(data = box_anno(comp2, var = "biomass", group_by = "RCP",
                            y = 600),
            aes(x, y, label = RCP, fill = NULL)) +
  labs(y = lab_bio0,
       x = NULL) +
  # code below so that legend doesn't appear but figure dimension are unchanged
  scale_fill_manual(values = cols_pft5_other,
                    guide = guide_legend(override.aes = list(alpha = 0))) +
  theme(legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),
        axis.text.x = element_blank()) 

# % of total biomass (stacked)
g2 <- g +
  geom_bar(aes(y = bio_perc), position = 'stack', stat = 'identity') +
  geom_text(data = box_anno(comp2, var = "bio_perc", group_by = "RCP",
                            y = 108),
            aes(x, y, label = RCP, fill = NULL)) +
  labs(y = lab_perc_bio0) +
  scale_fill_manual(values = cols_pft5_other, name = "") +
  scale_y_continuous(breaks = seq(0, 100, by = 25))


jpeg("figures/biomass/comp_stacked-bar_RCP8.5-M_c4on.jpeg",
     res = 600, width = 4, height = 6, units = 'in')
grid::grid.draw(
  # doing this so the left axis is aligned, and figures are the same height
  rbind(ggplotGrob(g1), ggplotGrob(g2))
)
dev.off()

# fire --------------------------------------------------------------------

# boxplot of fire return interval by RCP, and grazing intensity

map(levs_c4, function(lev_c4) {
  g1 <- fire1 %>% 
    filter(c4 == lev_c4) %>% 
    group_by(id) %>% 
    # note--this method of filtering before drawing boxplot 
    # creates biased boxplots
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
    scale_x_discrete(labels = id2year) +
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
    # NOTE--should update code to use 
    # compute_boxplot_stats() and stat = "identity" instead of remove_outliers()
    remove_outliers(var = "fire_return_diff") %>% 
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
    scale_x_discrete(labels = id2year) +
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

