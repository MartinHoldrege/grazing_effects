# Martin Holdrege

# script started Jan. 7, 2022

# Purpose--various functions used when creating maps
# This script sourced by other functions


# dependencies ------------------------------------------------------------

source("src/fig_params.R") # for cols_map_bio() color ramp function
source("src/fig_functions.R")

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

# for calculating the 2nd lowest and 2nd highest values
# across 13 GCMs
calc_low <- function(x) {
  
  if("SpatRaster" %in% class(x)) {
    stopifnot(nlyr(x) == 13)
    out <- terra::sort(x, decreasing = FALSE)[[2]]
    names(out) <- 'low'
  } else if (is.numeric(x)) {
    
    stopifnot(length(x) == 13 | length(x) == 1)
    
    if(length(x) == 13) {
      out <- sort(x, decreasing = FALSE, na.last = TRUE)[[2]]
    } else {
      out <- NA
    }
  } else {
    stop('class not supported')
  }

  out
}

calc_high <- function(x) {
  
  if("SpatRaster" %in% class(x)) {
    stopifnot(nlyr(x) == 13)
    out <- terra::sort(x, decreasing = TRUE)[[2]]
    names(out) <- 'low'
  } else if (is.numeric(x)) {
    
    stopifnot(length(x) == 13 | length(x) == 1)
    
    if(length(x) == 13) {
      out <- sort(x, decreasing = TRUE, na.last = TRUE)[[2]]
    } else {
      out <- NA
    }
    
  } else {
    stop('class not supported')
  }
  
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

#' get range from raster
#'
#' @param x spatraster with one or more layers
#' @param absolute get absolute range (for difference maps)
range_raster <- function(x, absolute = FALSE) {
  out <- terra::minmax(x) %>% 
    unlist() %>% 
    range(na.rm = TRUE)
  
  if(absolute) {
    m <- max(abs(out))
    out <- c(-m, m)
  }
  out
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
plot_map2 <- function(r, add_coords = TRUE,...)  {
  
  
  s <- stars::st_as_stars(r)
  
  map <- newRR3::plot_map(s, 
                  st_geom_state = states,
                  add_coords = add_coords,
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
                           values = NULL,
                           midpoint = NULL,
                           oob = scales::squish
)  {
  
  
  limits_inset <- if(is.null(limits))  {
    c(NA, NA)
  }  else 
    limits
  
  inset <- newRR3::inset_densitycountplot(as.numeric(values(r)),
                                  limits = limits_inset,
                                  add_vertical0 = add_vertical0) +
    scale_x_continuous(n.breaks = 4)
  
  s <- stars::st_as_stars(r)
  
  # allows for asymetrical colors (i.e. if midpoint isn't in the middle)
  if(!is.null(midpoint) & is.null(values)) {
    values <- values_about_midpoint(colors = colors, midpoint = midpoint,
                                    limits = limits)
    
  }
  
  map <- newRR3::plot_map(s, 
                  st_geom_state = states,
                  add_coords = TRUE) +
    newRR3::ggplot2_map_theme() +
    scale_fill_gradientn(na.value = 'transparent',
                         limits = limits,
                         name = scale_name,
                         colors = colors,
                         values = values,
                         oob = oob) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), units = 'in')) +
    newRR3::add_tag_as_label(tag_label) 
  
  map + inset_element2(inset)
  
}



#' 20 panel map
#'
#' @param r raster with 20 layers
#' @param info dataframe from create_rast_info
#' @param type_absolute value in the 'type' column of info that denotes layers of absolute
#' (non change) values
#' @param type_diff  value in the type column denoting change layers
#' @param title figure title
#' @param name4absolute variable name for the 'absolute' layers
#' @param name4diff variable name for the 'difference' layers
#' @param legend_title_absolute title for the absolute layers
#' @param legend_title_diff title for the difference layers
#' @param palette_absolute color pallete
#' @param palette_diff color palette
#'
#' @returns
#' 4x5 plot matrix
plot_map_20panel <- function(
    r,
    info,
    type_absolute,
    type_diff,
    title = NULL,
    name4absolute = "",
    name4diff = NULL,
    legend_title_absolute = name4absolute,
    legend_title_diff = name4diff,
    palette_absolute = cols_map_bio(10),
    palette_diff = cols_map_bio_d2,
    lims_diff = NULL,
    midpoint_diff = NULL
) {
  
  if(is.null(name4diff)) {
    name4diff <- paste("\u0394", name4absolute)
  }
  if(is.null(legend_title_absolute)) {
    legend_title_absolute <-  name4absolute
  }
  
  if(is.null(legend_title_diff)) {
    legend_title_diff <-  name4diff
  }
  
  
  # setup lookups
  types <- c(type_absolute, type_diff)
  names4leftside <- c(name4absolute, name4diff)
  names(names4leftside) <- types
  names4legend <- list(legend_title_absolute,legend_title_diff)
  names(names4legend) <- types
  cols <- list(palette_absolute, palette_diff)
  names(cols) <- types
  midpoints <- list(NULL, midpoint_diff)
  names(midpoints) <- types
  
  stopifnot(
    c(types) %in% info$type)
  
  info <- info %>% 
    mutate(type = factor(.data$type, levels = c(!!type_absolute, !!type_diff))) %>% 
    arrange(type, RCP, years, graze) %>% 
    mutate(tag_label = fig_letters[1:length(type)])
  
  
  stopifnot(nrow(info) == 4*5,
            is.factor(info$graze),
            unique(info$graze) == levels(info$graze),
            length(levels(info$graze)) == 4,
            info$id %in% names(r))
  
  graze_labs <- paste0(levels(info$graze), '\nGrazing')
  
  label_left <- info %>% 
    filter(.data$graze == levels(.data$graze)[1]) %>% 
    mutate(rcp_year = rcp_label(.data$RCP, .data$years,
                                include_parenth = FALSE),
           label_left = paste0(rcp_year, '\n', names4leftside[type])) %>% 
    pull(label_left)
  
  stopifnot(length(label_left) == 5)
  
  lims_absolute <- range_raster(r[[info$id[info$type == type_absolute]]], 
                       absolute = FALSE)
  
  if(is.null(lims_diff)) {
    lims_diff <- range_raster(r[[info$id[info$type == type_diff]]], absolute = TRUE)
  }
  lims <- list(lims_absolute, lims_diff)
  names(lims) <- types
  
  plots1 <- pmap(info[c('id', 'type', 'tag_label')], 
                 function(id, type, tag_label) {
                   plot_map_inset(r = r[[id]],
                                  colors = cols[[type]],
                                  tag_label = tag_label,
                                  limits = lims[[type]],
                                  scale_name = names4legend[[type]],
                                  midpoint = midpoints[[type]])
                 })
  
  # labels for the top margin
  plots_top <- map(graze_labs, function(x) {
    ggplot() +
      theme_void() +
      annotate("text", x = 1, y = 1, label = x, size = 3) +
      theme(plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
  })
  
  # labels for the left margin
  plots_left <- map(label_left, function(x) {
    ggplot() +
      theme_void() +
      annotate("text", x = 1, y = 1, label = x, angle = 90, size = 3) +
      theme(plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
  })
  
  plots2 <- c(list(plot_spacer()), plots_top)
  
  start <- 1
  for (i in 1:length(plots_left)) {
    end <- start + 3
    plots2 <- c(plots2, plots_left[i], plots1[start:end])
    start <- end + 1
  }
  
  g <- patchwork::wrap_plots(plots2, nrow = 6,
                             widths = c(0.18, rep(1,4)), 
                             heights = c(0.18, rep(1, 5))) + 
    plot_layout(guides = 'collect')
  
  g2 <- g&theme(legend.position = 'bottom')
  
  if(!is.null(title)) {
    g2 <- g2 + plot_annotation(subtitle = title)
  }
  
  g2
}


# crs -----------------------------------------------------------------------

# the crs to be used for the sagebrush conservation design (this is the same
# one as used by NLCD)
crs_scd <- terra::crs("PROJCRS[\"Albers_Conical_Equal_Area\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Albers Equal Area\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",23,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-96,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")



# polygons for basemaps ---------------------------------------------------

states <- sf::st_as_sf(spData::us_states) %>% 
  sf::st_transform(crs = crs_scd)

load_wafwa_ecoregions <- function() {
  shp1 <- sf::st_read("../SEI/data_raw/files_from_DaveT/WAFWAecoregionsFinal.shp")
  
  shp2 <- sf::st_transform(shp1, crs = crs_scd)
  
  shp2$ecoregion <- c("Great Basin", "Intermountain", "Plains")
  shp2
}
