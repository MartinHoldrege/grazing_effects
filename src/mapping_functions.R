# Martin Holdrege

# script started Jan. 7, 2022

# Purpose--various functions used when creating maps
# This script sourced by other functions


# dependencies ------------------------------------------------------------

source("src/fig_params.R") # for cols_map_bio() color ramp function

# breaks ------------------------------------------------------------------

#' create breaks and color vectors for plotting
#'
#' @param x Numeric vector which are the cell values from the raster
#' @param n_breaks Number of breaks to create
#' @param cols Vector of colors to use
#' @param n_legend_breaks number of breaks to show in the legend
#'
#' @return A list that contains `legendbks` which are the breaks to show
#' in the legend. `truebks` which are the breaks to pass the plotting function
#' (the have n_break number of breaks)
#' `cols` which are the colors to use for the legend, and `truecols` which
#' are the colors to use for the plotting function. `zero`, is a logical
#' of whether to minimum value was 0.
#'
#' @examples
#' create_breaks_cols(1:20)
#' create_breaks_cols(-1:20)
#' create_breaks_cols(0:20)
create_breaks_cols <- function(x, 
                          n_breaks = 17, 
                          cols = cols_map_bio(n_breaks - 1),
                          n_legend_breaks = 9
                          ) {
  stopifnot(
    is.numeric(x),
    n_breaks - 1 == length(cols) # should be 1 more break than colors
  )
  
  x <- x[!is.na(x)]
  
  # rounding so both the plotting function and legend use break 
  # points at whole numbers
  min <- floor(min(x))
  
  max <- ceiling(max(x))
  
  zero <- FALSE
  if(min == 0) {
    min <- floor(min(x[x != 0]))
    zero <- TRUE
  }
  
  # breaks
  bks <- seq(from = min, to = max, length.out = n_breaks)
  
  # breaks for the legend
  legendbks <- seq(from = min, to = max, length.out = n_legend_breaks)
  
  if (zero) {
    # making vector of where plotting function should break at, subtracting
    # small amount from 2nd element, so that 'min' doesn't get binned
    # into the grey zero category. 
    truebks <- c(-0.1, bks[1] - 0.001, bks[-1]) 
    truecols <- c("darkgrey", cols) # adding separate category for zeros
    # if zero is the minimum value (id for biomass, not for delta biomass)
  } else {
    truebks <- bks
    truecols <- cols
  }
  
  out <- list(truebks = truebks,
              legendbks= legendbks,
              cols = cols, 
              truecols = truecols,
              zero = zero)
  
  out
}


# operations on rasters ---------------------------------------------------

#' Scaled percent change, between two rasters
#' 
#' @description Calculates the percent change  of a target
#' layer relative to a reference layer (scaled by the max). 
#'
#' @param rast A SpatRast (terra) object with multiple layers
#' @param ref_layer Character (name) or number of the reference layer,
#' if length > 1, it should be the same length as the target_layer vector
#' @param target_layer Character (name) or number of the target layer
#'
#' @return SpatRast object (with as many layers as the length of the
#' target_layer vector)
rast_diff <- function(rast, ref_layer, target_layer) {
  
  # maximum value in the reference layer(s)
  max <- global(rast[[ref_layer]], max, na.rm = TRUE)$max
  
  if(nlyr(rast[[ref_layer]]) > 1 & 
     nlyr(rast[[ref_layer]]) != nlyr(rast[[target_layer]])) {
    stop("ref_layer needs to have same length as target layer, or have length of 1")
  }
  
   # scaled percent change
  out <- (rast[[target_layer]] - rast[[ref_layer]])/max*100 
  out
}


# maps --------------------------------------------------------------------

#' Create map of biomass data in western states
#' 
#' @description Note: this was designed to use SpatRaster objects, so 
#' the 'terra' package should be loaded before using the function
#'
#' @param rast Raster of class "SpatRaster", from the terra package. 
#' @param subset Number or string that identifies the layer of the rast
#' that should be plotted
#' @param title Title of figure
#' @param vec Numeric vector (optional), to base the breaks/cut points used
#' in plotting on. By default it uses the values of the cells in the raster.
#' @param show0legend logical, if false suppresses the appearance of the 
#' legend that shows 0 values in grey
#' @param legend_lab label for the main legend
#' @param n_breaks number of color bins used
#'
#' @return A map of biomass
image_bio <- function(rast, subset, title = "", vec = NULL, 
                      show0legend = TRUE,
                      legend_lab = expression("Biomass ("*gm^-2*")"),
                      n_breaks = 17 
                      ) {
  
  stopifnot(
    length(subset) == 1 # this function can only work with one raster layer
  )
  
  # extract raster cell data as a fector
  if(is.null(vec)) {
    # I suspect extracting the actual values into memory like
    # this is inefficient
    vec <- values(subset(rast, subset = subset))
  }
  
  # create list of colors and break points to use
  b <- create_breaks_cols(x = vec, n_breaks = n_breaks)
  
  # main figure
  image(subset(rast, subset = subset), # the layer to be plotted
        # if default used, maps are more pixelated looking, not that with this
        # higher resolution an error is sometimes thrown the first time the funcion
        # is called (a similar issue
        # is documented here: https://github.com/rspatial/terra/issues/30)
        maxcell = 500000,
        col = b$truecols, 
        breaks = b$truebks, 
        ylim = c(30, 49),
        xlim = c(-125, -102.7), useRaster = TRUE,
        xlab = "", ylab ="",
        bty = "n", xaxt = "n",yaxt="n")
  mtext(title, side = 3, line = 0, adj = 0, cex=0.7)
  maps::map("state", interior = T, add = T)
  
  # Legend for zero values
  if(show0legend & b$zero) {
    polygon(x = c(-125,-125,-119,-119), y = c(32,35,35,32),
            col = "white",border = "white")
    legend(x = -120, y = 35.5, xjust = 0.7,
           legend = expression("Biomass ("*gm^-2*") = 0"), 
           fill = "darkgrey", bty = "n",cex = 0.8, 
           border = "darkgrey", x.intersp = 0.2)
  }
  
  # Color bar/legend at the bottom
  polygon(x = c(-125,-125,-102.7,-102.7), y = c(29,32.7,32.7,29), 
          border = "white", col = "white")
  rasterImage(t(b$cols), -125, 31.75, -102.7, 32.5)
  polygon(x = c(-125, -125, -102.7 ,-102.7), y = c(31.75, 32.5, 32.5, 31.75), 
          lwd = 1.5)
  mtext(legend_lab, side = 1, line = -0.49, cex = 0.7)
  axis(side = 1, pos = 31.75, at = seq(-125,-102.7, length.out = 9),
       cex.axis = 0.9, labels = b$legendbks)
  
}

#' Create map of % change of biomass in western states
#' 
#' @description At the moment this figure uses fixed cut points for the colors
#'
#' @param rast Raster of class "SpatRaster", from the terra package. 
#' @param subset Number or string that identifies the layer of the rast
#' that should be plotted
#' @param title Title of figure
#'
#' @return Map of percent change in biomass
image_bio_diff <- function(rast, subset, title = "") {
  
  # range of bio_diff is: -85.67349 103.51214
  stopifnot(
    length(subset) == 1 # this function can only work with one raster layer
  )
  

  # for now hard coding breaks and colors
  truebks <- sort(c(-(2^(1:7)), 0, 2^(1:7)))
  cols <-  cols_map_bio_d
  
  # run a check on whether truebks are excluding values
  vec <- values(subset(rast, subset = subset))
  min <- min(vec, na.rm = TRUE)
  max <- max(vec, na.rm = TRUE)
  if(min < min(truebks) | max > max(truebks)) {
    stop("Raster contains data outside of the range of the breaks")
  }

  
  # main figure
  image(subset(rast, subset = subset), # the layer to be plotted
        # if default used, maps are more pixelated looking
        maxcell = 500000, 
        col = cols, 
        breaks = truebks, 
        ylim = c(30, 49),
        xlim = c(-125, -102.7), useRaster = T,
        xlab = "", ylab ="",
        bty = "n", xaxt = "n",yaxt="n")
  mtext(title, side = 3, line = 0, adj = 0, cex= 0.7)
  maps::map("state", interior = T, add = T)
  
  # Color bar/legend at the bottom
  polygon(x = c(-117,-117,-109,-109), y = c(32,34.5,34.5,32),col = "white",border = "white")
  polygon(x = c(-125,-125,-102.7,-102.7), y = c(30,33.5,33.5,30),border = "white", col = "white")
  
  # creating mini polygons for each legend color
  axisat <- seq(-125, -102.7, length.out = length(cols) + 1)
  for (i in 1:length(cols)){
    polygon(x = c(axisat[i],axisat[i],axisat[i+1],axisat[i+1]), 
            y = c(31.75,32.5,32.5,31.75),border = cols[i], col = cols[i])
    
  }
  # black outline of color bar
  polygon(x = c(-125, -125, -102.7 ,-102.7), y = c(31.75, 32.5, 32.5, 31.75), 
          lwd = 1.5)
  mtext(expression(paste("", Delta, " Biomass (%)")), side = 1, line = -0.49, cex = 0.7)
  axis(side = 1, pos = 31.75, at = seq(-125,-102.7, length.out = 15),
       cex.axis = 0.9, labels = truebks)
  
}


#' Map of minimum grazing level needed to cross threshold
#'
#' @param rast raster
#' @param subset layer of raster to plot (integer or character string)
#' @param title map title (string)
#'
#' @return Map
image_min_gr <- function(rast, 
                         subset, 
                         title = ""
) {
  
  stopifnot(
    length(subset) == 1 # this function can only work with one raster layer
  )
  
  # 
  color <- c("darkgrey", rev(RColorBrewer::brewer.pal(4, "YlOrRd")))

  # main figure
  image(subset(rast, subset = subset), # the layer to be plotted
        # if default used, maps are more pixelated looking
        maxcell = 500000,
        col = color,
        breaks = seq(1:6) -0.5, 
        ylim = c(30, 49),
        xlim = c(-125, -102.7), useRaster = TRUE,
        xlab = "", ylab ="",
        bty = "n", xaxt = "n",yaxt="n")
  mtext(title, side = 3, line = 0, adj = 0, cex=0.7)
  maps::map("state", interior = T, add = T)
  
  legend(
    -125, 33.5,
    legend = c("Light grazing", 
               "Moderate grazing", 
               "Heavy grazing", 
               "Very Heavy Grazing",
               "Didn't cross threshold"),
    title = "Grazing required to cross threshold",
    fill = color,
    bg = "white",
    box.col = "white",
    ncol = 2,
    cex = 0.7)
}


