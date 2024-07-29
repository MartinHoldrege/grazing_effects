# Martin Holdrege

# Script started 1/6/2022

# Purpose--to create figures of biomass data from simulations
# This script does not create maps of the interpolated data. 

# Note--many of these boxplots have repeated code, some of that could be 
# consolidated into functions

# dependencies ------------------------------------------------------------

library(tidyverse)
library(lemon) # for facet_rep_wrap function
library(patchwork)
source("scripts/02_summarize_bio.R")
source("src/fig_params.R")
source("src/fig_functions.R") # box_ann function defined here

# params ------------------------------------------------------------------

theme_set(theme_custom1())
line_loc2 <- c(1.5, 3.5) # for figures w/ 2 vertical lines

# width of figs
wfig_box1 <- 9 # width of boxplots in inches
wfig_box2 <- 6 # width of boxplots that show just 5 pft groups
hfig_box2 <- 5 # height of 5 panel boxplots

# number of columns of panels in boxplots
ncol_box <- 3

# * vectors/dfs for 'looping' ---------------------------------------------

# for looping
levs_pft <- unique(pft5_bio2$PFT) # all pft levels in the main datafile

# for when need to 'loop' over both pft4 and c4 on and off
levs_pft_run <- expand_grid(pft = factor(levs_pft, levs_pft), 
                           run = runs_graze) %>% 
  arrange(pft, run)

levs_grefs_run <- expand_grid(ref_graze = names(pft5_d_grefs), 
                             run = runs_graze)
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

# boxplots of utilization and biomass for all functional groups
pdf("figures/biomass/bio-util-boxplot_all.pdf",
    height = 8, width = wfig_box1)

figs_l <- map(runs_graze, function(x) {
  out <- list()
  df <-  filter(pft5_bio_b, run == x)
  
  out$bio <- ggplot(df, aes(x = id, y = biomass, fill = RCP)) +
    box1() +
    labs(caption = x)
  
  out$util <- ggplot(df, aes(x = id, y = utilization, fill = RCP)) +
    box1(var = 'utilization',
         y = lab_util0) +
    labs(caption = x)
  out
})
map(figs_l, \(x) x$bio)
map(figs_l, \(x) x$util)
dev.off()

# c4 on and off in side by side panels 
# for each pft
# code removed [see pre '24 commit]

# ** fewer scenarios ------------------------------------------------------

# showing the biomass boxplot but only for the main PFTs, and RCP mid-century
# for possible use in manuscript

df <- pft5_bio_b %>% 
  filter_rcp_run(PFT = TRUE)

n_zero <- df %>% 
  group_by(RCP, PFT, graze) %>% 
  summarize(n_zero = sum(biomass == 0 | is.na(biomass), na.rm =TRUE),
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

g0 <- ggplot(df, aes(x = RCP, fill = graze)) +
  facet_rep_wrap(~ PFT, scales = "free", ncol = 3) +
  scale_fill_graze() +
  labs(x = lab_rcp) +
  theme_box_pft5()

suffix <- filter_rcp_run(df = NULL)

jpeg(paste0("figures/biomass/pub_qual/bio-box_pft5_", suffix, ".jpeg"),
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")
g0 +
  labs(y = lab_bio0) +
  geom_text(data = n_zero2, aes(y = y, label = string),
            size = 2,vjust = 'inward',
            lineheight = 0.7) +  
  geom_boxplot(aes(y = biomass), outlier.size = outlier.size) +
  geom_vline(xintercept = 1.5, linetype = 2)

dev.off()

jpeg(paste0("figures/biomass/pub_qual/util-box_pft5_", suffix, ".jpeg"),
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")

g0 +
  labs(y = lab_util0) +
  geom_boxplot(aes(y = utilization), outlier.size = outlier.size) +
  geom_vline(xintercept = 1.5, linetype = 2)

dev.off()

# * scatterplot (vs climate) ----------------------------------------------

# each combination PFt and C4 on or off is plotted on seperate pages
# facets are RCP/year combinations, colors are grazing
pdf("figures/biomass/bio_vs_climate_v2.pdf",
    width = 6, height = 5)

# here pmap is working by taking the column names (i.e. iterating over rows)
pmap(levs_pft_run, function(pft, run) {
  pft <- as.character(pft)

  df <- pft5_bio2[pft5_bio2$PFT == pft & pft5_bio2$run == run, ] 
  
  base <- function(ylab){
    list(
      facet_rep_wrap(~RCP + years),
                 scale_color_graze(),
      labs(y = ylab,
           subtitle = pft,
           caption = paste(run, '\n x-axis shows current climate in all panels')),
        theme(legend.position = c(0.85, 0.15),
              axis.text = element_text(size = 7)) 
      )
    }
  
  g1 <- ggplot(df, aes(y = biomass, color = graze)) +
    base(ylab = lab_bio0)

  g2 <- ggplot(df, aes(y = utilization, color = graze)) +
    base(ylab = lab_util0)
  
  list(bio = climate_scatter(g1), util = climate_scatter(g2))
  
})

# facets are grazing level, colors are RCP
pmap(levs_pft_run, function(pft, run) {
  pft <- as.character(pft)
  df <- pft5_bio2[pft5_bio2$PFT == pft & pft5_bio2$run == run &
                       pft5_bio2$years != "2030-2060", ] 
  
  base <- function(ylab){
    list(
      facet_rep_wrap(~graze),
      scale_color_manual(values = cols_rcp),
      labs(y = ylab,
           subtitle = pft,
           caption = paste(run,
                           "\nOnly showing 2070-2100 for RCP4.5 & 8.5")),
      theme(legend.position = "top",
              axis.text = element_text(size = 7)) 
    )
  }
  
  g1 <- ggplot(df, aes(y = biomass, color = RCP)) +
    base(ylab = lab_bio0)
  
  g2 <- ggplot(df, aes(y = utilization, color = RCP)) +
    base(ylab = lab_util0)
  
  list(bio = climate_scatter(g1), util = climate_scatter(g2))

})

dev.off()

# * cheat light vs heavy graze --------------------------------------------
# comparing cheatgrass light vs heavy grazing, to help determine why some
# sites have such large % changes in cheatgrass (it appears it is the
# group of sites with low cheat biomass)

cheat_df <- pft5_bio2 %>% 
  filter(PFT == "Cheatgrass", RCP %in% c("Current", "RCP45")) %>% 
  select(run, years, RCP, graze, site, biomass, utilization)

pdf("figures/biomass/light_vs_heavy_graze_scatterplot.pdf",
    width = 6, height = 6)

for(r in runs_graze) {
df <- cheat_df %>% 
  filter(graze == "Light") %>% 
  left_join(filter(cheat_df, graze == "Heavy"),
            by = c("run", "RCP", "years", "site"),
            suffix = c("_light", "_heavy")) %>% 
  filter(run == r)

g <- ggplot(df) +
  facet_rep_wrap(~RCP + years, ncol = 2) +
  geom_abline(slope = 1) +
  labs(caption = r)

print(g +
  labs(x = "Biomass (light grazing)",
       y = "Biomass (heavy grazing)") +
  geom_point(aes(biomass_light, biomass_heavy)))

print(g +
  labs(x = "Utilization (light grazing)",
       y = "Utilization (heavy grazing)") +
  geom_point(aes(utilization_light, utilization_heavy))
  )
}
dev.off()
# compare runs  -----------------------------------------------------------

# * absolute biomass ------------------------------------------------------

# select which two runs to compare with each other
rn1 <- runs_graze[1]
rn2 <- runs_graze[2]

pft5_wide <- pft5_bio2 %>% 
  filter(run %in% runs_graze) %>% 
  pivot_wider(
    id_cols = c('years', 'RCP', 'graze', 'id', 'PFT', 'site'),
    names_from = 'run',
    values_from = c('biomass', 'utilization')
  )

compare_runs_scatter <- function(pft, var) {
  pft5_wide %>% 
    filter(.data$PFT == pft) %>% 
    ggplot(aes(.data[[paste0(var, '_', rn1)]], 
               .data[[paste0(var, '_', rn2)]], 
               color = graze)) +
    geom_point(alpha = 0.6, size = 1) +
    facet_rep_wrap(RCP ~ years) +
    geom_abline(slope = 1, color = "black") +
    geom_smooth(method = "lm") +
    labs(
      x = rn1,
      y = rn2,
      caption = paste(
        "Regression lines seperate for each grazing treatment,",
        "1:1 line in black."),
      subtitle = paste(pft, var, "at the site level, comparing runs")
    ) +
    scale_color_graze() +
    theme(legend.position = c(0.85, 0.15))
}

pdf("figures/biomass/bio-util_compare-runs.pdf",
    width = 9, height = 7)

map(levs_pft, function(x) {
  list(compare_runs_scatter(x, var = 'biomass'), 
       compare_runs_scatter(x, var = 'utilization'))
})

dev.off()

# C3Pgrass/Pgrass ---------------------------------------------------------
g <- C3_Pgrass_ratio %>% 
  filter(run %in% runs_graze) %>% 
  ggplot(aes(id, C3_Pgrass_ratio, fill = RCP)) +
  box1(add_facet_wrap = FALSE, var = "C3_Pgrass_ratio",
       group_by = c("run", "graze"),
       y = lab_C3_Pgrass_ratio) +
  facet_rep_wrap(~run) +
  labs(subtitle = "Ratio of C3Pgrass biomass to total Pgrass biomass",
       caption = "Separately showing data from simulations with C4 expansion on and off") 
g
ggsave("figures/biomass/C3_Pgrass_ratio.jpeg", g, 
       width = 6, height = 4)

# biomass change -------------------------------------------------------
# response variable is change in biomass

# * change relative to same graze ------------------------------------

# boxplot of change in biomass (scaled percent and effect size), 
# for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period, relative to current
# conditions for that same grazing level (i.e. this is the climate effect)


pdf("figures/biomass/pft5_bio-util_diff_boxplots_v1.pdf",
    height = 8, width = wfig_box1)

# scaled % change (see pre-2024 code)

# effect size
l <- map(runs_graze, function(rn) {
  
  df <- pft5_bio_es1 %>% 
    filter(run == rn)
  
  g1 <- ggplot(df, aes(id2, bio_es, fill = graze)) +
    box2(axis_data = pft5_bio_es1, 
         var = "bio_es", repeat.tick.labels = "y") +
    add_sec_axis() +
    labs(y = lab_es0,
         caption = rn)
  
  g2 <- ggplot(df, aes(id2, util_es, fill = graze)) +
    box2(axis_data = pft5_bio_es1, 
         var = "bio_es", repeat.tick.labels = "y") +
    add_sec_axis() +
    labs(y = lab_es0_util,
         caption = rn)
  
  list(bio = g1, util = g2)
})

l

dev.off()

# * change relative to reference graze ------------------------------------


# ** boxplot ---------------------------------------------------------------

# boxplots showing change in biomass (and effect size) relative to current time 
# period and a given grazing intensity ('reference class')

pdf("figures/biomass/pft5_bio-util_diff_gref_boxplots_v1.pdf", 
    height = 6.5, width = wfig_box1)

# % change
# see pre 2024 commits for % change

# effect size
pmap(levs_grefs_run, function(ref_graze, run){
  rn <- run
  df0 <- pft5_es_grefs[[ref_graze]]
  
  df <- df0 %>% 
    filter(run == rn) 
  
  g1 <- ggplot(df, aes(id2, bio_es, fill = graze)) +
    box3(axis_data = df0, var = "bio_es") +
    add_sec_axis() +
    labs(y = lab_es0,
         subtitle = paste("Change in biomass relative to", tolower(ref_graze), 
                          "grazing \n under current conditions"),
         caption = rn)
  
  g2 <- ggplot(df, aes(id2, util_es, fill = graze)) +
    box3(axis_data = df0, var = "util_es") +
    add_sec_axis() +
    labs(y = lab_es0_util,
         subtitle = paste("Change in utilization relative to", tolower(ref_graze), 
                          "grazing \n under current conditions"),
         caption = rn)
  list(g1, g2)
})


dev.off()


# ** scatterplot (vs climate) ----------------------------------------------

# effect size vs MAP and MAT, for for ref class of light grazing, for
# each RCP/time period and PFT
pdf("figures/biomass/bio-util-diff_vs_climate_v1.pdf",
    width = 6, height = 5)
# scaled %change[see pre 2024 code]

# effect size
# just showing results for default
l2 <- pmap(levs_pft_run[levs_pft_run$run == runs_graze['default'], ], 
           function(pft, run) {
  # using Light grazing as reference class
  df <- pft5_es_grefs[["Light"]] %>% 
    filter(PFT == as.character(pft))
  
  g1 <- df %>% 
    filter(.data$run == !!run) %>% 
    ggplot(aes(y = bio_es, color = graze)) +
    scatter_light(pft, run, axis_data = df,
                  subtitle_response = 'biomass') +
    labs(y = lab_es0) +
    add_sec_axis()
  
  g2 <- df %>%
    filter(.data$run == !!run) %>%
    ggplot(aes(y = util_es, color = graze)) +
    scatter_light(pft, run, axis_data = df,
                  subtitle_response = 'utilization') +
    labs(y = lab_es0_util) +
    add_sec_axis()
  
  if(pft == 'Sagebrush') {
    # utilization is 0 for sagebrush for light grazing, so figure doesn't work
    out <- list(climate_scatter(g1))
  } else {
    out <- list(climate_scatter(g1),
                climate_scatter(g2))
  }
  out
})

l2

dev.off()


# * wgcm ----------------------------------------------

# change in biomass going from light grazing, to some other grazing level,
# for a given scenario, ie. change within a gcm/climate scenario
# (e.g. light grazing RCP 8.5 end century to heavy graze, RCP 8.5 end century)


# ** boxplots -------------------------------------------------------------
pdf("figures/biomass/pft5_bio_diff_wgcm_boxplots_v2.pdf", 
    height = 8, width = wfig_box1)

# effect size
map(runs_graze, function(x){
 
  out <- list()
  # w/outliers
  df <- pft5_es_wgcm %>% 
    filter(run == x) 
  out[[1]] <- ggplot(df, aes(id2, bio_es, fill = graze)) +
    box2(axis_data = pft5_es_wgcm,
         var = "bio_es", xintercept = line_loc2,
         subtitle = "Reference group is light grazing for the given climate scenario",
         repeat.tick.labels = TRUE) +
    add_sec_axis() +
    labs(y = lab_es0,
         caption = paste0("All data shown\n", x)) +
    theme(axis.text.y.right = element_text(size = 5)) 

   
   # throws warnings b/ 0 utilization for sagebrush light grazing
   out[[2]] <- ggplot(df, aes(id2, util_es, fill = graze)) +
     box2(axis_data = pft5_es_wgcm, 
          var = "util_es", xintercept = line_loc2,
          subtitle = "Reference group is light grazing for the given climate scenario",
          repeat.tick.labels = TRUE) +
     add_sec_axis() +
     labs(y = lab_es0_util,
          caption = paste0("All data shown\n", x)) +
     theme(axis.text.y.right = element_text(size = 5)) 
  
  # see pre-'24 commits for version of  boxplot w/ outliers removed
  out
})

dev.off()


# ***fewer scenarios ------------------------------------------------------

# effect size boxplot/violin plot for main PFTs, and only RCP8.5 mid-century,
# (i.e. closer to publication quality)

df <- pft5_es_wgcm_heavy %>% 
  filter_rcp_run(PFT = TRUE)


# base of the plot
g <- ggplot(df, aes(x = RCP, fill = graze)) +
  geom_blank() + # added so that geom_vline doesn't throw an error
  geom_vline(xintercept = 1.5, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.3) +
  facet_rep_wrap(~ PFT, ncol = 3,
                 repeat.tick.labels = FALSE) +
  scale_fill_graze(exclude = 'Heavy') +
  labs(x = lab_rcp) +
  theme_box_pft5()+
  add_sec_axis(name = "% Change relative to heavy grazing") +
  scale_color_manual(values = cols_graze) +
  guides(color = 'none')


suffix <- filter_rcp_run(df = NULL)
jpeg(paste0("figures/biomass/pub_qual/bio-diff_boxplot_pft5_", suffix, ".jpeg"),
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")
g+
  geom_boxplot(aes(y = bio_es),outlier.size = outlier.size) +
  theme(panel.spacing.x = unit(-2, "lines"))+
  labs(y = lab_es1_heavy)

dev.off()

jpeg(paste0("figures/biomass/pub_qual/util-diff_boxplot_pft5_", suffix, ".jpeg"),
     res = 600, height = hfig_box2, width = wfig_box2, units = "in")
g+
  geom_boxplot(aes(y = util_es),outlier.size = outlier.size) +
  theme(panel.spacing.x = unit(-2, "lines"))+
  labs(y = lab_es1_util_heavy)
  
dev.off()


# % below threshold ------------------------------------------------------

ref_threshold2 <- ref_threshold %>% 
  mutate(string = paste(round(threshold, 0), "g/m^2"))

pdf("figures/threshold/threshold_dotplots_v2.pdf", 
    height = 6.5, width = wfig_box1)

map(runs_graze, function(x) {
  
threshold1 %>% 
  filter(run == x) %>% 
  ggplot(aes(x = id, y = pcent_above)) +
  geom_text(data = ~box_anno(., var = "pcent_above",
                             group_by = c("PFT", "graze"),
                             y = 110),
            aes(x, y, label = graze, fill = NULL),
            size = 2.5) +
  # printing the value of the threshold on each panel
  geom_text(data = filter(ref_threshold2, run == x),
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
       caption = paste0(x,
                       "\nValue of threshold provided in each panel"))
  
})

dev.off()

# * fewer scenarios -------------------------------------------------------
# threshold dotplot just showing current and RCP 8.5 mid-century.
# grouped by climate scenario.

df <- threshold1 %>% 
  filter_rcp_run(PFT = TRUE) %>% 
  # it's not relevant whether cheatgrass drops below a 'sustainable' level.
  filter(PFT != "Cheatgrass" ) 

ref_threshold3 <- ref_threshold2 %>% 
  filter(run %in% df$run, PFT %in% df$PFT)

jpeg(paste0("figures/threshold/threshold_5th-pcent_", suffix, ".jpeg"),
     res = 600, units = 'in',
     height = 4, width = 4)
ggplot(df, aes(x = RCP, y = pcent_above)) +
  facet_rep_wrap(~ PFT, ncol = 2) +
  geom_point(aes(color = graze), position = position_dodge(width = 1)) +
  scale_color_graze() +
  geom_vline(xintercept = 1.5, linetype = 2) +
  labs(x = lab_rcp,
       y = lab_below0) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8)) +
  geom_text(data = ref_threshold3,
            aes(x = 0.8, y = 25, label = string),
            size = 2.5)
dev.off()

# composition -------------------------------------------------------------

# stacked barcharts showing the composition of the 5 main PFTs.
# this is based on the mean biomass across sites for each of the PFTs.
# also separate barcharts for only the herbaceous PFTs

# function for creating the stacked bar charts below
g_comp <- function(df, color, text_y = 600, ylab2 = lab_perc_bio0,
                   text_y2 = 108, # height of text for composition figure
                   shift_ylim = 0,
                   group_by = 'RCP', label = "RCP",
                   xintercept = 4.5) {
  g <- ggplot(df, aes(x = id, fill = PFT)) +
    scale_x_discrete(labels = id2graze) +
    labs(x = lab_graze) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) +
    geom_vline(xintercept = xintercept, linetype = 2) 
  
  # absolute biomass (stacked)
  g1 <- g +
    geom_bar(aes(y = biomass), position = 'stack', stat = 'identity') +
    geom_text(data = box_anno(df, var = "biomass", group_by = group_by,
                              y = text_y),
              aes_string("x", "y", label = label, fill = NULL)) +
    labs(y = lab_bio0,
         x = NULL) +
    # code below so that legend doesn't appear but figure dimension are unchanged
    scale_fill_manual(values = color,
                      guide = guide_legend(override.aes = list(alpha = 0))) +
    theme(legend.title = element_text(color = "transparent"),
          legend.text = element_text(color = "transparent"),
          axis.text.x = element_blank()) +
    coord_cartesian(ylim = c(0, text_y + text_y*shift_ylim))
  
  # % of total biomass (stacked)
  g2 <- g +
    geom_bar(aes(y = bio_perc), position = 'stack', stat = 'identity') +
    geom_text(data = box_anno(df, var = "bio_perc", group_by = group_by,
                              y = text_y2),
              aes_string("x", "y", label = label, fill = NULL)) +
    labs(y = ylab2) +
    scale_fill_manual(values = color, name = "") +
    scale_y_continuous(breaks = seq(0, 100, by = 25),
                       limits = c(0, text_y2 + text_y2*shift_ylim)) 
    
  
  # doing this so the left axis is aligned, and figures are the same height
  rbind(ggplotGrob(g1), ggplotGrob(g2))
}

jpeg(paste0("figures/biomass/comp_stacked-bar_", suffix, ".jpeg"),
     width = 8, height = 6, units = 'in', res = 600)

grid::grid.draw(
  cbind(g_comp(comp2, color =  cols_pft5_other),# All pfts
        g_comp(comp_herb1, color =  cols_herb, text_y = 180,
               ylab2 = "% of biomass"))
)

dev.off()


# * all scenarios ---------------------------------------------------------
# stacked bar charts for the appendix that show all climate scenarios. 

# all plant functional types
comp2_all2 <- comp2_all %>% 
  ungroup() %>% 
  arrange(RCP, years, graze) %>% 
  mutate(
    id = paste(RCP, years, graze, sep = "_"),
    id = factor(id, levels = unique(id)),
    years = ifelse(years == 'Current', "", as.character(years)),
         RCP_years = paste(RCP, years, sep = "\n"))

g <- g_comp(comp2_all2, color =  cols_pft5_other, label = "RCP_years", 
       group_by = "RCP_years", xintercept = 0.5 + 4*1:4, text_y = 630,
       text_y2 = 113, shift_ylim = 0.05)

jpeg(paste0("figures/biomass/comp_stacked-bar_pft5_all-scenarios_", 
            runs_graze['default'], ".jpeg"),
     width = 8, height = 6, units = 'in', res = 600)

grid::grid.draw(g)

dev.off()

# herbaceous functional type
# composition of just key herbaceous PFTs
comp_herb_all <- comp2_all2 %>% 
  filter(PFT %in% c('C3Pgrass', 'C4Pgrass', 'Pforb', "Cheatgrass")) %>% 
  # this df is already grouped
  group_by(id) %>% 
  mutate(bio_perc = biomass/(sum(biomass))*100)

g2 <- g_comp(comp_herb_all, color =  cols_pft5_other, label = "RCP_years", 
             group_by = "RCP_years", xintercept = 0.5 + 4*1:4, text_y = 190,
             text_y2 = 113, shift_ylim = 0.05)

jpeg(paste0("figures/biomass/comp_stacked-bar_herb_all-scenarios_", 
            runs_graze['default'], ".jpeg"),
     width = 8, height = 6, units = 'in', res = 600)

grid::grid.draw(g2)

dev.off()

# fire --------------------------------------------------------------------

# boxplot of fire return interval by RCP, and grazing intensity

map(runs_graze, function(x) {
  df <- fire_med1 %>% 
    filter(run == x) 
  g1 <- ggplot(df, aes(x = id, y = fire_prob, fill = RCP)) +
    geom_text(data = ~box_anno(., var = "fire_prob",
                               group_by = c("graze")),
              aes(x, y, label = graze, fill = NULL),
              size = 2.5) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_rcp, name = "Scenario") +
    scale_x_discrete(labels = id2year) +
    geom_vline(xintercept = line_loc, linetype = 2) +
    theme(legend.position = "top",
          axis.text = element_text(size = 7))  +
    labs(x = lab_yrs,
         y = lab_firep0) 
  
  g2 <- ggplot(df, aes(x = rcp_label(RCP, years, include_parenth = FALSE), 
                 y = fire_prob, fill = graze))+
    geom_boxplot() +
    scale_fill_graze() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
    labs(y = lab_firep0,
         x = NULL)
  
  # fire difference boxplot
  g3 <- fire_d_wgraze %>% 
    filter(run == x) %>% 
    group_by(id) %>% 
    #remove_outliers(var = "fire_prob_diff") %>% 
    ggplot(aes(id2, fire_prob_diff, fill = graze)) +
    geom_text(data = ~box_anno(., var = "fire_prob_diff", 
                               group_by = c("RCP"), id = "id2"),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_boxplot(position = "dodge",
                 outlier.color = NA) +
    scale_fill_graze() +
    # so just display the RCP
    scale_x_discrete(labels = id2year) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 2.5, linetype = 2) +
    labs(x = lab_yrs,
         y = lab_firep1)
  
  jpeg(paste0("figures/fire/fire_return_boxplots_", x, ".jpeg"),
       height = 7, width = 7, res = 600, units = "in")
  print(g1 + g2 + g3 + guide_area() + patchwork::plot_layout(guides = 'collect'))
  dev.off()
})

