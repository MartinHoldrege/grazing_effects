# Martin Holdrege

# Script started: March 8, 2022

# Purpose: 
# Create ternary plots of the composition of 3 plant functional types at 
# a time, and composition changes under climate/grazing treatments


# dependencies ------------------------------------------------------------

library(ggtern)
source("scripts/02_summarize_bio.R")
source("src/fig_params.R")
source("src/fig_functions.R")
theme_set(theme_classic())


# functions ---------------------------------------------------------------


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
create_lines_df <- function(df1, df2, xyz) {
  
  stopifnot(nrow(df1) == 200,
            nrow(df2) == 200,
            xyz %in% names(df1),
            xyz %in% names(df2))
  
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
create_lines_df2 <- function(df, split_var, xyz) {
  df_list <- split(df, df[[split_var]], drop = TRUE)
  stopifnot(length(df_list) == 2)
  df1 <- df_list[[1]]
  df2 <- df_list[[2]]
  
  create_lines_df(df1, df2, xyz = xyz)
}


# * fig functions ---------------------------------------------------------


# function to create terniary plots
tern <- function(df, split_var, xyz) {
  ggtern(df, aes_string(x = xyz[1], y = xyz[2], z = xyz[3])) +
    geom_segment(data = create_lines_df2(df, split_var, xyz = xyz), 
                 aes(x, y, z, xend = xend, yend = yend, zend = zend),
                 alpha = 0.2) +
    theme_arrowlarge() +
    labs(x = "", y = "", z="",
         xarrow = paste(xyz[1], "cover"),
         yarrow = paste(xyz[2], "cover"),
         zarrow = paste(xyz[3], "cover"))
}

# function to create ternary plots that compare light vs heavy grazing
# in a given climate scenario
tern_wrcp <- function(df, 
                      xyz, 
                      levs_rcp = levs_2rcp, 
                      c4string = 'c4on') {
  
  out <- map(levs_rcp, function(x){
    # c4on and off are the same for current conditions here,
    # so not making a fiture in that case
    if (x == "Current" & c4string == 'c4off') return (NULL)
    
    df <- df %>% 
      filter(RCP == x,
             c4 == c4string) 
    
    lab <- if (x!='Current') {
      paste0(x, ' (', unique(df$years), ')')
    } else {
      x
    }
    g <- tern(df, split_var = 'graze', xyz = xyz) +
      geom_point(aes(color = graze), size = 0.5) +
      scale_color_manual(values = cols_graze[levs_2graze]) +
      labs(subtitle = paste0("Light vs heavy grazing, under ",
                             lab, " conditions (", c4string, ")"))
    g
  })
  out
}

# ternary plots show current vs RCP8.5, under given grazing level
tern_wgraze <- function(df, # dataframe
                        xyz, # string--names of 3 PFTs of interest
                        levs_graze = levs_2graze, # grazing levels of interest
                        c4string = 'c4on') {
  df <- df %>% 
    filter(c4 == c4string)
  
  # creating ternary plots for each level of two grazing levels
  out <- map(levs_graze, function(x) {
    df %>% 
      filter(graze == x) %>% 
      tern(split_var = 'RCP', xyz = xyz) +
      geom_point(aes(color = RCP), size = 0.5) +
      scale_color_manual(values = cols_rcp[levs_2rcp]) +
      labs(subtitle = paste0("Current vs RCP8.5 (2030-2060), under ",
                             tolower(x), " grazing (", c4string, ")"))
  })
  out
}

# create ternary plots comparing current light grazing 
# vs future heavy grazing
tern_across_graze_RCP <- function(df, xyz, c4string = 'c4on') {
  df %>% 
    filter((RCP == "RCP8.5" & graze == 'Heavy')|
             (RCP == "Current" & graze == "Light"),
           c4 == c4string)%>% 
    tern(split_var = 'graze', xyz = xyz) +
    geom_point(aes(color = graze, shape = RCP), size = 0.5) +
    scale_color_manual(values = cols_graze[levs_2graze]) +
    labs(subtitle = paste0("Current light vs RCP8.5 (2030-2060) heavy grazing",
                           " (", c4string, ")"))
}

# preparing data ----------------------------------------------------------

pft_levs <- list(
  Pherb = c("Pforb", "C3Pgrass", "C4Pgrass"), # herbaceous perennials
  SC34 = c("Sagebrush", "C3Pgrass", "C4Pgrass"), # sagebrush and c3/c4 grasses
  SGF = c("Sagebrush", "Pgrass", "Pforb") # sage and grass and forbs
)

# create dataframes for each set of 3 PFTs, that contain % of biomass
# (where the percent of biomass of the 3 pfts sums to 100)
perc_dfs <- map(pft_levs, function(x) {
  out <- pft5_bio2 %>% 
    select(-id) %>% 
    filter(PFT %in% x,
           RCP %in% c("Current", "RCP8.5"),
           years %in% c("Current", "2030-2060"),
           graze %in% c("Light", "Heavy")) %>% 
    group_by(c4, years, RCP, graze, site) %>% 
    mutate(perc = biomass/sum(biomass)*100) %>% 
    select(-biomass, -indivs) %>% 
    pivot_wider(names_from = "PFT", values_from = "perc")
  out
})


# figures -----------------------------------------------------------------


pdf("figures/ternary/ternary_by-PFT-group_v1.pdf")

# * within RCP level comparisons --------------------------------------------
map2(perc_dfs, pft_levs, function(df, xyz) {
  out <- list()
  out$c4on <- tern_wrcp(df, xyz = xyz, c4string = 'c4on')
  out$c4off <- tern_wrcp(df, xyz = xyz, c4string = 'c4off')
  out
})

# * within grazing level comparisons ----------------------------------------

# ie ternary plots show current vs RCP8.5, under given grazing level

map2(perc_dfs, pft_levs, function(df, xyz) {
  out <- list()
  out$c4on <- tern_wgraze(df, xyz = xyz, c4string = 'c4on')
  out$c4off <- tern_wgraze(df, xyz = xyz, c4string = 'c4off')
  out
})

# * across gazing and RCP comparison ----------------------------------------

# comparing current light grazing 
# vs future heavy grazing

# looping over pft groupings
map2(perc_dfs, pft_levs, function(df, xyz) {
  out <- list()
  out$c4on <- tern_across_graze_RCP(df, xyz = xyz, c4string = 'c4on')
  out$c4off <- tern_across_graze_RCP(df, xyz = xyz, c4string = 'c4off')
  out
})


dev.off()


