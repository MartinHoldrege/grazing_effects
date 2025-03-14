# Martin Holdrege

# script started Jan. 7, 2022

# Purpose--various functions used when creating maps
# This script sourced by other functions


# dependencies ------------------------------------------------------------

source("src/fig_params.R") # for cols_map_bio() color ramp function

# misc --------------------------------------------------------------------

# the crs to be used for the sagebrush conservation design (this is the same
# one as used by NLCD)
crs_scd <- terra::crs("PROJCRS[\"Albers_Conical_Equal_Area\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Albers Equal Area\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",23,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-96,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")


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

#' fill a template with values
#' 
#' @description
#' Filling a spatraster instead of just using values() because
#' data (df) doesn't include rows for missing values. 
#' 
#'
#' @param template template raster (terra spatraster)
#' @param df dataframe with 1 or more data columns and 1 cellnum column
#' where the cellnums correspond to the cellnumbers of the template raster
#' @param use_template_values logical, use the actuall cell velues of the
#' template rasters as the cell numbers to match by (as opposed to just
#' calculating the cellnumbers of the template). Important distinction
#' if the template raster was created, and contains cell numbers but was
#' then trimmed so cell values don't match the cellnumbers derived from ncell
#' 
#' @return spatraster
fill_raster <- function(df, template,
                        use_template_values = FALSE) {
  stopifnot(is.data.frame(df),
            'cell_num' %in% names(df) | 'cellnumber' %in% names(df))
  if('cellnumber' %in% names(df)) {
    df <- rename(df, cell_num = cellnumber)
  }
  
  if(!use_template_values) {
    full_df <- tibble(cell_num = 1:terra::ncell(template[[1]]))
    
    full_df2 <- full_join(full_df, df, by = 'cell_num') %>% 
      dplyr::select(-cell_num)
    
    stopifnot(nrow(full_df) == nrow(full_df2))
  } else {
    stopifnot(terra::nlyr(terra::rast(template)) == 1)
    full_df <- as.data.frame(template) %>% 
      tidyr::drop_na()
    names(full_df) <- 'cell_num'
    
    stopifnot(all(full_df$cell_num %in% df$cell_num))

    full_df2 <- left_join(full_df, df, by = 'cell_num') %>% 
      dplyr::select(-cell_num)
  }

  
  r_out <- terra::rast(template[[1]], nlyrs = ncol(full_df2))
  names(r_out) <- names(full_df2)
  
  if(!use_template_values) {
    terra::values(r_out) <- as.matrix(full_df2)
  } else {
    r_out[!is.na(r_out)] <- as.matrix(full_df2)
  }
  
  r_out
}

# maps --------------------------------------------------------------------


inset_element2 <- function(x) {
  patchwork::inset_element(
    x,
    0.005, 0.005, 360 / 1133, 230 / 1236, # left, bottom, right, top in npc units
    align_to = "panel",
    clip = TRUE,
    ignore_tag = TRUE
  )
}

# some of these functions rely on the newRR3 package (by Daniel Schlaepfer)
plot_map2 <- function(r, ...)  {
  
  
  s <- stars::st_as_stars(r)
  
  map <- newRR3::plot_map(s, 
                  st_geom_state = states,
                  add_coords = TRUE,
                  ...) +
    newRR3::ggplot2_map_theme() +
    theme(plot.tag = element_text(face = 'plain'))
  
  map
  
}

plot_map_inset <- function(r,
                           colors = colors,
                           tag_label = "",
                           scale_name = NULL,
                           limits = NULL,
                           add_vertical0 = FALSE,
                           values = NULL
)  {
  
  
  limits_inset <- if(is.null(limits))  {
    c(NA, NA)
  }  else 
    limits
  
  inset <- newRR3::inset_densitycountplot(as.numeric(values(r)),
                                  limits = limits_inset,
                                  add_vertical0 = add_vertical0)
  
  s <- stars::st_as_stars(r)
  
  map <- newRR3::plot_map(s, 
                  st_geom_state = states,
                  add_coords = TRUE) +
    newRR3::ggplot2_map_theme() +
    scale_fill_gradientn(na.value = 'transparent',
                         limits = limits,
                         name = scale_name,
                         colors = colors,
                         values = values) +
    newRR3::add_tag_as_label(tag_label) 
  
  map + inset_element2(inset)
  
}

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

#' legend for bio diff map
#'
#' @param cols vector (colors to use)
#' @param truebks vector of where breaks will occur
#' @param cex defines size of legend label
#'
#' @return plots the legend on top of the current graphic
legend_bio_diff <- function(cols, truebks, cex = 0.7) {
  
  # making background white
  # polygon(x = c(-117,-117,-109,-109), y = c(32,34.5,34.5,32),col = "white",
  #         border = "white")
  polygon(x = c(-125,-125,-102.7,-102.7), y = c(-29.8,33.5,33.5,30),
          border = "white", col = "white")
  
  # Color bar/legend at the bottom
  # creating mini polygons for each legend color
  axisat <- seq(-125, -102.7, length.out = length(cols) + 1)
  for (i in 1:length(cols)){
    polygon(x = c(axisat[i],axisat[i],axisat[i+1],axisat[i+1]), 
            y = c(31.75,32.5,32.5,31.75),border = cols[i], col = cols[i])
    
  }
  # black outline of color bar
  polygon(x = c(-125, -125, -102.7 ,-102.7), y = c(31.75, 32.5, 32.5, 31.75), 
          lwd = 1.5)
  mtext(expression(paste("", Delta, " Biomass (%)")), side = 1, line = -0.49, 
        cex = cex)
  axis(side = 1, pos = 31.75, at = seq(-125,-102.7, length.out = length(truebks)),
       cex.axis = 0.9, labels = truebks)
}

#' Create map of % change of biomass in western states
#' 
#' @description At the moment this figure uses fixed cut points for the colors
#'
#' @param rast Raster of class "SpatRaster", from the terra package. 
#' @param subset Number or string that identifies the layer of the rast
#' that should be plotted
#' @param title Title of figure
#' @param legend logical, whether to add a legend to the plot
#' @param adjust_ylim logical, if true then adjust ylim when there is no
#' legend (i.e. if no legend the ylims are restricted)
#' @param cex defines size of legend label
#'
#' @return Map of percent change in biomass
image_bio_diff <- function(rast, subset, title = "", legend = TRUE,
                           adjust_ylim = TRUE, cex = 0.7) {
  
  # range of bio_diff is: -85.67349 103.51214
  stopifnot(
    length(subset) == 1, # this function can only work with one raster layer
    is.logical(legend)
  )
  
  # for now hard coding breaks and colors
  bks <- c(60, 40, 30, 20, 10, 5)
  # last break is 150 so that larger increases are included
  # consider improving
  truebks <- c(-100, -bks, 0, rev(bks), 150)
  cols <-  cols_map_bio_d
  
  # run a check on whether truebks are excluding values
  vec <- values(subset(rast, subset = subset))
  min <- min(vec, na.rm = TRUE)
  max <- max(vec, na.rm = TRUE)
  if(min < min(truebks) | max > max(truebks)) {
    warning(title, " Raster contains data outside of the range of the breaks",
         'max value is ', round(max), ". Min value is ", round(min), '.')
  }

  if (legend | !adjust_ylim) {
    ylim <- c(30, 49) # ylim larger to accommodate legend
  } else {
    ylim <- c(33.5, 49)
  }
  
  # main figure
  image(subset(rast, subset = subset), # the layer to be plotted
        # if default used, maps are more pixelated looking
        maxcell = 500000, 
        col = cols, 
        breaks = truebks, 
        ylim = ylim,
        xlim = c(-125, -102.7), useRaster = T,
        xlab = "", ylab ="",
        bty = "n", xaxt = "n",yaxt="n")
  mtext(title, side = 3, line = 0, adj = 0, cex= 0.7)
  maps::map("state", interior = T, add = T)
  
  if(!adjust_ylim) {
    # making background white
    # polygon(x = c(-117,-117,-109,-109), y = c(32,34.5,34.5,32),col = "white",
    #         border = "white")
    polygon(x = c(-125,-125,-102.7,-102.7), y = c(-29.8,33.5,33.5,30),
            border = "white", col = "white")
  }
  
  # add legend
  if (legend) {
    legend_bio_diff(cols = cols, truebks = truebks, cex = cex)
  }
  
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
  
  color <- rev(RColorBrewer::brewer.pal(9, "YlOrRd"))[c(1, 3, 5, 7)]
  color <- c(color, "#abd9e9")
  
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

# crs -----------------------------------------------------------------------

# the crs to be used for the sagebrush conservation design (this is the same
# one as used by NLCD)
crs_scd <- terra::crs("PROJCRS[\"Albers_Conical_Equal_Area\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Albers Equal Area\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",23,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-96,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")



# polygons for basemaps ---------------------------------------------------

states <- sf::st_as_sf(spData::us_states) %>% 
  sf::st_transform(crs = crs_scd)


