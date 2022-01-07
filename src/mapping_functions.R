# Martin Holdrege

# script started Jan. 7, 2022

# Purpose--various functions used when creating maps
# This script sourced by other functions

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


# maps --------------------------------------------------------------------

image_bio <- function(rast, subset, title = "", vec = NULL) {
  
  stopifnot(
    length(subset) == 1 # this function can only work with one raster layer
  )
  
  # extract raster cell data as a fector
  if(is.null(vec)) {
    vec <- values(subset(rast, subset = subset))
  }
  
  # create list of colors and break points to use
  b <- create_breaks_cols(x = vec, n_breaks = 17)
  
  # main figure
  image(rast, 
        y = subset, # the layer to be plotted
        col = b$truecols, 
        breaks = b$truebks, 
        ylim = c(30, 49),
        xlim = c(-125, -102.7), useRaster = T,
        xlab = "", ylab ="",
        bty = "n", xaxt = "n",yaxt="n")
  mtext(title, side = 3, line = 0, adj = 0, cex=1)
  maps::map("state", interior = T, add = T)
  
  # Legend for zero values
  if(b$zero) {
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
  mtext(lab_bio0, side = 1, line = -0.49, cex = 0.7)
  axis(side = 1, pos = 31.75, at = seq(-125,-102.7, length.out = 9),
       cex.axis = 0.9, labels = b$legendbks)
  
}

