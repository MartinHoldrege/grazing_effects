# Martin Holdrege

# Script started March 16, 2022

# Purpose: 
# Create maps showing the composition of 3 PFTs, by assigning
# proportion of total biomass to red, green, and blue. Percent biomass
# rasters values are scaled between 0, and 1 across related rasters
# so for example, use the min and max % biomass across all sagebrush
# grazing levels to scale, if comparing grazing levels.

# Note--currently this script can't be run as a job
# due to this error: Error in x@ptr$nrow() : external pointer is not valid

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(tmap)
library(Ternary)
library(spData) # for us_states polygon
source("src/general_functions.R")

# user defined values -----------------------------------------------------

# change depending on whether want to run this code for c4 expansion
# functionality turned on or off
# Note--this functionality isn't fully working at the moment for c4off
# (i think because c4off current wasn't upscaled (b/ current c4on and 
# off are the same))
c4string <- "c4on" # "c4off" # 

# read in data ------------------------------------------------------------

p <- "data_processed/interpolated_rasters"

# *rasters ----------------------------------------------------------------

# median biomass across GCMs
med2 <- rast(file.path(p, "bio_future_median_across_GCMs.tif"))

# * info about rasters ----------------------------------------------------

# names of the raster layers, and treatment info, 
rast_info_med <- readRDS(file.path(p, "raster_info.RDS")) %>% 
  select(-layer_num, -GCM, -id) %>% 
  filter_rcp_c4(c4string = c4string) %>% 
  filter(graze %in% c("Light", "Heavy")) %>% 
  distinct() 

# functions --------------------------------------------------

# move these functions to there own .R file, if they need
# to be used in another script

# calculate percent of total biomass, in each of 3 layers 
# of a SpatRast raster object
calc_perc <- function(x) {
  stopifnot(class(x) %in% c("raster", "SpatRaster"),
            nlyr(x) == 3 # must have 3 layers (later for r, g, b)
  )
  total <- sum(x) # total biomass across the 3 layers
  
  out <- x/total*100 # calculate percent of total
  out
}


#' Adjust min and max vectors
#' 
#' @description Goal is to to create min/max values so that each
#' vertex of the ternary plot sums to 100%. Currently this is 
#' achieved by adjust the max values (ie increasing them as needed)
#'
#' @param min vector of length 3
#' @param max vector of length 3
#'
#' @return list with two elements, min and max
#' 
adjust_minmax <- function(min, max) {
  
  stopifnot(length(min) == 3,
            length(max) == 3, # one value for each of 3 layers
            # should be a percent
            min <= 100 & min >= 0,
            max <= 100 & max >= 0,
            letters[1:3] %in% names(min),
            letters[1:3] %in% names(max))
  
  # top vertex
  max['a'] <- 100 - sum(min[c('b', 'c')])
  # bottom right vertex
  max['b'] <- 100 - sum(min[c('a', 'c')])
  # bottom left vertex
  max['c'] <- 100 - sum(min[c('a', 'b')])
  
  # check that each vertex sums to 100
  stopifnot(
    sum(max['a'], min[c("b", "c")]) == 100,
    sum(max['b'], min[c("a", "c")]) == 100,
    sum(max['c'], min[c("a", "b")]) == 100
  )
  
  list(min = min, max = max)
}

#' Normalize SpatRaster
#'
#' @param x SpatRaster object
#' @param min vector, min values to use for each layer of x, when normalizing
#' @param max vector, max values to use for each layer of x, when normalizing
#'
#' @return SpatRaster object where values range from 0 to 1
rast_normalize <- function(x, min = NULL, max = NULL) {
  stopifnot(class(x) == "SpatRaster")
  
  n <- nlyr(x)
  
  # if vector not provided calculate from raster
  if(is.null(min)) {
    min <- t(minmax(x))[, 1] # min values for each layer
    min <- floor(min)
  } else {
    stopifnot(length(min)==n)
  }
  
  if(is.null(max)) {
    min <- t(minmax(x))[, 2] # min values for each layer
    max <- floor(max)
  } else {
    stopifnot(length(max)==n)
  }
  
  out <- (x - min)/(max - min)
  
  # check that normalization worked
  minmax_out <- unname(t(minmax(out)))
  stopifnot(minmax_out[, 1] >= 0, # min should be 0 (or greater)
            minmax_out[, 2] <= 1 # max should be 1 (or less)
  )
  
  out
}

#' Add grid lines to ternary plot
#'
#' @description Draws 4 lines on ternary plots in each direction
#' This if for drawing on top of colors (at which point the underlying
#' colors get overplotted and aren't visible)
add_lines <- function() {
  breaks = c(0.2, 0.4, 0.6, 0.8) # where lines will be created
  breaks2 <- rep(breaks, each = 2)
  
  rows <- c(1, 4, 5, 8)
  # horizontal lines (a-b lines)
  m <- matrix(0, nrow = length(breaks)*2, ncol = 3)
  m[ ,1] <- rep(rev(breaks), each = 2)
  
  m[-rows, 2] <- breaks2[-rows]
  m[rows, 3] <- breaks2[rows]
  AddToTernary(lines, m, col = 'darkgrey', lty = 'dotted', lwd = 2)
  # vertical a to c lines
  m2 <- m[nrow(m):1, c(3, 2, 1)]
  AddToTernary(lines, m2, col = 'darkgrey', lty = 'dotted', lwd = 2)
  # vertical b to c lines
  m3 <- m[, c(2, 1, 3)]
  AddToTernary(lines, m3, col = 'darkgrey', lty = 'dotted', lwd = 2)
  
  # draw triangle
  m4 <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1,0, 0), nrow = 3)
  AddToTernary(lines, m4, col = 'darkgrey', lty = 'solid', lwd = 3)
}

# create rgb triangle
ternary_rgb <- function(alab, blab, clab, # labels for three axis
                        # min and max axis labels for the 3 labels:
                        amin, amax,
                        bmin, bmax,
                        cmin, cmax,
                        shift = 0.05 # how much to move labels outward
) {
  TernaryPlot(alab = paste(alab, "\u2192"), 
              blab = paste(blab, "\u2192"),
              clab = paste("\u2190", clab),
              lab.col = c('red', 'darkgreen', 'blue'),
              #point = 'right', 
              lab.cex = 0.8,
              grid.minor.lines = 0,
              grid.lines = 5,
              grid.lty = 'solid',
              col = rgb(1, 1, 1),
              grid.col = 'darkgrey',
              axis.col = rgb(0.6, 0.6, 0.6),
              ticks.col = rgb(0.6, 0.6, 0.6),
              axis.rotate = FALSE,
              axis.labels = FALSE,
              padding = 0.08
  )
  
  cols <- TernaryPointValues(rgb)
  ColourTernary(cols, spectrum = NULL)
  add_lines() # add grid lines
  
  cex = 0.6
  
  # add labels for min/max values of the axes
  text(-0.5 - shift, 0, amin, cex = cex)
  text(0 - shift, sqrt(1-0.5^2), amax, cex = cex)
  text(0 + shift, sqrt(1-0.5^2) + shift, bmin, cex = cex)
  text(0.5+ shift, 0 + shift, bmax, cex = cex)
  text(0.5 + shift, 0-shift, cmin, cex = cex)
  text(-0.5 + shift*0.5, 0 - shift, cmax, cex = cex)
  
}

# test
ternary_rgb("a", "b", "c",
            "amin", "amax",
            "bmin", "bmax", 
            "cmin", "cmax",
            shift = 0.06)

# create a temporary png file of rgb triangle,
# to be read in and added as a 'logo' when making
# a map with tmap. Note this function returns the path to the
# temp file created
temp_ternary_rgb <- function(pfts, # character vector length 3
                             min, # numeric vector length 3
                             max # numeric vector lenght 3
) {
  # path for creating a temporary file
  # including pfts in name, to reduce errors later
  s <- paste(pfts, collapse = "_") 
  temp <- paste0(tempfile(s), ".png")
  
  abc <- c("a", "b", "c")
  
  # vectors need to be named
  stopifnot(abc %in% names(pfts),
            abc %in% names(min),
            abc %in% names(max))
  
  pfts <- paste("%", pfts) # for labeling
  names(pfts) <- abc # names get removed when pasting
  
  # save temp file 
  png(temp, height = 2, width = 2, units = "in",
      res = 300)
  par(oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0))
  
  # creating rgb triangle
  ternary_rgb(pfts['a'],
              pfts['b'],
              pfts['c'],
              amin = min['a'],
              amax = max['a'],
              bmin = min['b'],
              bmax = max['b'],
              cmin = min['c'],
              cmax = max['c']
  )
  dev.off()
  
  temp # return path to the temp png file
}

# function to create rgb transition maps.
# designed to be passed to map()
create_rgb_map <- function(RCP, # string
                          graze, # string
                          data, # dataframe, with 3 rows (one for each layer of interest in
                          # the raster)
                          # it must include columns labeled id_noGCM, years, and c4
                          rast_perc, # raster (with 3 layers, which gives the percent of
                          # the biomass, represented by the 3 PFTS
                          min, # minimum value to use for normalization
                          max, # max value to use for normalization
                          temp # path to the temporary rgb triangle png)
) {
  
  # normalize (so values fall between 0 and 1)
  r <- rast_perc[[data$id_noGCM]]
  r <- rast_normalize(r, min = min, max = max)
  
  # test that temp file is at least for the correct pfts
  s <- paste(data$PFT, collapse = "_")
  if(!str_detect(temp, s)){
    stop("temp file doesn't match the desired PFTs")
  }
  
  # create title title
  if(RCP != 'Current') {
    # adding in time period
    RCP <- paste0(RCP, "(", data$years, ")")
  }
  
  c4 <- data$c4 %>% unique()
  title <- paste0(RCP, ", ", graze, " grazing (", c4, ")")
  
  # create map
  out <- tm_logo(temp,
                 position = c("LEFT", "BOTTOM"), height = 6,
                 margin = 0) +
    tm_shape(r) +
    tm_rgb(max.value = 1,
           legend.show = TRUE) +
    tm_shape(us_states) +
    tm_borders() +
    tm_layout(main.title = title,
              main.title.size = 1)
  out
}



# maps --------------------------------------------------------------------


# * example map -----------------------------------------------------------

pfts <- c(a = "Pforb", b = "C3Pgrass", c = "C4Pgrass")

id <- paste0("c4on_", pfts, "_biomass_Current_Current_Light")

x <- calc_perc(med2[[id]])

# later on the goal is to use min and max derived from a larger
# group of rasters, so multiple figures have the same color triangles.

min <- t(minmax(x))[, 1]
max <- t(minmax(x))[, 2]
names(min) <- letters[1:3]
names(max) <- letters[1:3]


min <- floor(min)
max <- ceiling(max)

# adjust so vertices will sum to 100%
min <- adjust_minmax(min, max)$min # at the moment this step doesn't make any changes
max <- adjust_minmax(min, max)$max

r <- rast_normalize(x, min = min, max = max)

temp <- temp_ternary_rgb(pfts, min = min, max = max)

tm_logo(temp,
        position = c("LEFT", "BOTTOM"), height = 7,
        margin = 0) +
  tm_shape(r) +
  tm_rgb(max.value = 1,
         legend.show = TRUE) +
  tm_shape(us_states) +
  tm_borders() +
  tm_layout(main.title = "test title")


# * maps by pft group---------------------------------------------------------

pft_levs <- list(
  # herbaceous perennials
  Pherb = c(a = "Pforb", b = "C3Pgrass", c = "C4Pgrass"),
  # sagebrush and c3/c4 grasses
  SC34 = c(a = "Sagebrush", b = "C3Pgrass", c = "C4Pgrass"), 
  # sage and grass and forbs
  SGF = c(a = "Sagebrush", b = "Pgrass", c = "Pforb") 
)

# iterating over pft grouping,
# and creating 4 maps for each grouping
figs_by_pft <- map(pft_levs, function(pfts) {
  
  print(pfts)
  # info about the raster layers of interest
  info1 <- rast_info_med %>% 
    filter(PFT %in% pfts) %>% 
    # so ordering of PFTs is consistent
    mutate(PFT = factor(PFT, levels = pfts)) %>% 
    arrange(RCP, graze, PFT) 
  
  info2 <- info1 %>% 
    group_by(RCP, graze) %>% 
    nest()
  
  rast_perc <- map(info2$data, function(x) {
    r <- med2[[x$id_noGCM]] # 3 raster layers
    calc_perc(r) # calc percent of total for each of the 3 layers
  })
  rast_perc <- rast(rast_perc) # convert list back to SpatRast
  
  # the layer ids for each pft
  pft_id <- map(pfts, function(x) {
    id <- info1 %>% 
      filter(PFT == x) %>% 
      pull(id_noGCM)
  })
  
  # min for each pft, across rasters of interest for that pft
  pft_min <- map_dbl(pft_id, function(x) {
    out <- t(minmax(rast_perc[[x]]))[, 1] %>% 
      min() %>% 
      floor()
    out
  })
  
  # max 
  pft_max <- map_dbl(pft_id, function(x) {
    out <- t(minmax(rast_perc[[x]]))[, 2] %>% 
      max() %>% 
      ceiling()
    out
  })
  
  # adjust so vertices will sum to 100%
    # at the moment this function doesn't actually make changes to min
  pft_min <- adjust_minmax(pft_min, pft_max)$min 
  pft_max <- adjust_minmax(pft_min, pft_max)$max
  
  # rgb triangle
  temp <- temp_ternary_rgb(pfts, min = pft_min, max = pft_max)
  
  # create the maps (1 for each grazing/RCP scenario)
  out_figs <- pmap(info2, create_rgb_map,
                rast_perc = rast_perc, min = pft_min, max = pft_max,
                temp = temp)
  
  out_figs
})

# arrange figures into 2*x panels, 1 for each pft grouping
figs_arranged <- map(figs_by_pft, tmap_arrange, nrow = 2)

# save figures
file <- paste0("figures/composition_maps/composition_rgb_v2_", 
               c4string, ".pdf")

pdf(file, height = 9, width = 9)
  figs_arranged
dev.off()


