# Martin Holdrege

# Script started: March 8, 2022

# Purpose: 
# Create ternary plots of the composition of 3 plant functional types at 
# a time, and composition changes under climate/grazing treatments


# dependencies ------------------------------------------------------------

library(ggtern)
source("scripts/02_summarize_bio.R")
source("src/fig_params.R")
theme_set(theme_classic())


# functions ---------------------------------------------------------------

herb_levs <- c("C3Pgrass", "C4Pgrass", "Pforb")
levs_2graze <- c("Light", "Heavy")
levs_2rcp <- c("Current", "RCP8.5")


#' create dataframe for geom_segment 
#'
#' @param df1 dataframe of reference group values
#' @param df2 dataframe of comparison group values
#' @param xyz character vector of the 3 variables of interest 
#' (i.e. cover columns for 3 PFTs)
#'
#' @return combined dataframe based on df1, and df2, by site number
create_lines_df <- function(df1, df2, xyz = herb_levs) {
  
  stopifnot(nrow(df1) == 200,
            nrow(df2) == 200)
  df1a <- df1[, c('site', xyz)]
  names(df1a) <- c('site', 'x', 'y', 'z')
  df2a <- df2[, c('site', xyz)]
  names(df2a) <- c('site', 'xend', 'yend', 'zend')
  
  out <- full_join(df1a, df2a, by = 'site')
  stopifnot(nrow(out) == 200) # make sure join was good
  
  out
}

#' wrapper around create_lines_df
#'
#' @param df a dataframe
#' @param split_var column name to split by
#' @param xyz variables to plot
#'
#' @return dataframe ready to use in geom_segment() with ggtern()
create_lines_df2 <- function(df, split_var, xyz = herb_levs) {
  df_list <- split(df, df[[split_var]], drop = TRUE)
  stopifnot(length(df_list) == 2)
  df1 <- df_list[[1]]
  df2 <- df_list[[2]]
  
  create_lines_df(df1, df2, xyz = xyz)
}

# function to create terniary plots
tern <- function(df, split_var) {
  ggtern(df, aes(x = C3Pgrass, y = C4Pgrass, z = Pforb)) +
    geom_segment(data = create_lines_df2(df, split_var), 
                 aes(x, y, z, xend = xend, yend = yend, zend = zend),
                 alpha = 0.2) +
    theme_arrowlarge() +
    labs(x = "", y = "", z="",
         xarrow = "C3Pgrass cover",
         yarrow = 'C4Pgrass cover',
         zarrow = 'Pforb cover')
}

# preparing data ----------------------------------------------------------


# * grasses and forbs -----------------------------------------------------

# percent of total perennial herbaceous biomass belonging to each of the
# 3 perennial herbaceous PFTs
herb_perc1 <- pft5_bio2 %>% 
  filter(PFT %in% herb_levs,
         RCP %in% c("Current", "RCP8.5"),
         years %in% c("Current", "2030-2060"),
         graze %in% c("Light", "Heavy"),
         c4 == 'c4on') %>% 
  group_by(c4, years, RCP, graze, id, site) %>% 

  mutate(perc = biomass/sum(biomass)*100) %>% 
  select(-biomass, -indivs) %>% 
  pivot_wider(names_from = "PFT", values_from = "perc")

herb_perc1


# figures -----------------------------------------------------------------


pdf("figures/ternary/ternary_C3Pgrass-C4Pgass-Pforb_c4on_v1.pdf")

# * within RCP level comparisons --------------------------------------------

for (x in levs_2rcp) {
  df <- herb_perc1 %>% 
    filter(RCP == x) 
  lab <- if (x!='Current') {
    paste0(x, ' (', unique(df$years), ')')
  } else {
    x
  }
  g <- tern(df, split_var = 'graze') +
    geom_point(aes(color = graze), size = 0.5) +
    scale_color_manual(values = cols_graze[levs_2graze]) +
    labs(subtitle = paste("Light vs heavy grazing, under",
                          lab, "conditions"))
  print(g)
}


# * within grazing level comparisons ----------------------------------------


for (x in levs_2graze) {
  g <- herb_perc1 %>% 
    filter(graze == x) %>% 
    tern(split_var = 'RCP') +
    geom_point(aes(color = RCP), size = 0.5) +
    scale_color_manual(values = cols_rcp[levs_2rcp]) +
    labs(subtitle = paste("Current vs RCP8.5 mid-century, under",
      tolower(x), "grazing"))
  print(g)
}


# * across gazing and RCP comparison ----------------------------------------

herb_perc1 %>% 
  filter((RCP == "RCP8.5" & graze == 'Heavy')|
           (RCP == "Current" & graze == "Light"))%>% 
  tern(split_var = 'graze') +
  geom_point(aes(color = graze, shape = RCP), size = 0.5) +
  scale_color_manual(values = cols_graze[levs_2graze]) +
  labs(subtitle = "Current light vs RCP8.5 mid-century heavy grazing")

dev.off()


