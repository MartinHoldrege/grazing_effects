# Martin Holdrege

# script started Jan. 7, 2022

# Purpose--various functions used when creating maps
# This script sourced by other functions


# dependencies ------------------------------------------------------------

source("src/fig_params.R") # for cols_map_bio() color ramp function
source("src/fig_functions.R")
source("src/general_functions.R")
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
rast_diff_scaled <- function(rast, ref_layer, target_layer) {
  
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


#' calculate difference between historical and future scenarios
#'
#' @param r raster with  layers for current and future scenarios
#' @param info dataframe from create_rast_info()
#' @param type_absolute string--type part of the band names
#' @param type_diff new string to use to describe the 'difference' type
#' @param match_vars variables to match by when matching historical and future
#' layers
#' @param include_percent also included layers of % change in the output raster
rast_diff <- function(r, info, type_absolute, type_diff, 
         match_vars = c('run', 'PFT', 'type', 'graze', 'summary'), 
         include_percent = FALSE) {
  
  stopifnot(
    'Current' %in% info$RCP,
    length(unique(info$type)) == 1
  )
  
  info_cur <- info[info$RCP == 'Current', ]
  
  stopifnot(
    # duplicates mean a matching variable is likely missing
    sum(duplicated(info_cur[match_vars])) == 0
  )
  
  info_fut <- info[info$RCP != 'Current', ]
  
  info_comb <- left_join(info_fut, info_cur, by = match_vars,
                         suffix = c("_fut", "_cur"))
  
  diff <- r[[info_comb$id_fut]] - r[[info_comb$id_cur]]
  names(diff) <- stringr::str_replace(names(diff), type_absolute, type_diff)
  
  if(include_percent) {
    perc_diff <- diff/r[[info_comb$id_cur]]*100
    names(perc_diff) <- stringr::str_replace(names(perc_diff), 
                                             type_diff, 
                                             paste0(type_diff, '-perc'))
    out <- c(diff, perc_diff)
  } else {
    out <- diff
  }
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

#' calculate median (and low/high) summaries across GCMs
#'
#' @param r spatraster with GCM level results
#' @param info create_rast_info() output pdf, including id_noGCM column
#' @param include_low_high logical, also include low/high estimates
summarize_gcms_raster <- function(r, info, include_low_high) {
  
  stopifnot(c('id_noGCM', 'id') %in% names(info),
            names(r) %in% info$id,
            info$id %in% names(r))
  
  info <- arrange(info, .data$id)
  r <- r[[info$id]]# sort (makes later check easier)
  
  # list of rasters, each element is a SpatRaster including all GCMs
  # for a given treatment/scenario combination
  rast_gcm_l <- split(r, f = info$id_noGCM)
  
  names(rast_gcm_l) <- unique(info$id_noGCM)

  # check that the number of layers in each dataset is correct
  stopifnot(map_dbl(rast_gcm_l, nlyr) == info %>%
              group_by(id_noGCM) %>%
              summarise(n = n()) %>%
              pull(n)
  )
  
  # note these median rasters now inlcude the median under current
  # conditions (i.e. medians of one value). 
  med1 <- map(rast_gcm_l, app, fun = "median")
  
  med2 <- rast(med1)
  names(med2) <- paste0(names(med2), '_median')
  
  if(include_low_high) {
    # exclude current time-period
    rast_gcm_l_fut <- discard(rast_gcm_l, \(x) nlyr(x) == 1)
    # calculate low (2nd lowest) and high (2nd highest)
    # estimates across GCMs
    low <- map(rast_gcm_l_fut, calc_low)
    names(low) <- paste0(names(low), "_low")
    
    high <- map(rast_gcm_l_fut, calc_high)
    names(high) <- paste0(names(high), "_high")
    comb1 <- c(med2, rast(low), rast(high))
    comb1 <- comb1[[sort(names(comb1))]]
  } else {
    comb1 <- med2
  }
  
  comb1
}

#' calculate difference from current conditions
#'
#' @param r spatraster
#' @param info info dataframe from create_rast_info()
#' @param type_from regex for portion of file name to change (e.g. '_biomass')
#' @param type_to what to change it to (e.g. 'bio-rdiff-cref')
#' @param by vector of column names to join by (i.e. these are for matching
#' the right layers from current vs future)
calc_rast_cref <- function(r, info, type_from, type_to = NULL,
                            by = c('run', 'PFT', 'type', 'graze')) {
  
  if(is.null(type_to)) {
    type_to <- paste0(type_from, "-rdiff-cref")
  }
  
  stopifnot(
    c('id', 'RCP', by) %in% names(info),
    info$id %in% names(r)
  )
  info_cur <- info %>% 
    filter(.data$RCP == 'Current') 
  
  info_fut <- info %>% 
    filter(.data$RCP != 'Current')
  
  info_comb <- info_fut[c('id', by)] %>% 
    left_join(info_cur, suffix = c("_fut", "_cur"), by = by)
  
  diff <- r[[info_comb$id_fut]] - r[[info_comb$id_cur]]
  names(diff) <- str_replace(names(diff), type_from, type_to)
  diff
}

downsample <- function(r, n = 100) {
  spatSample(r, size = n, method = 'regular',as.raster = TRUE)
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
    theme(plot.tag = element_text(face = 'plain')) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
  
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
    values_absolute = NULL,
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
  values <- list(values_absolute, NULL)
  names(values)
  
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
                                  values = values[[type]],
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


# calculated percent change for 20panel maps (or any other)
calc_rast_perc_change <- function(r, info, type_absolute, type_diff) {
  info <- info %>% 
    mutate(type = factor(.data$type, levels = c(!!type_absolute, !!type_diff))) %>% 
    arrange(type, RCP, years, graze) 
  
  # make sure don't have multipe PFTs etc in this info/r
  stopifnot(sum(duplicated(info[c('type', 'RCP', 'years', 'graze')])) == 0)
  
  r_current <- r[[info$id[info$type == type_absolute & info$RCP == 'Current']]]
  r_future_diff <- r[[info$id[info$type == type_diff]]]
  r_perc <- r_future_diff/r_current*100
  c(r[[info$id[info$type == type_absolute]]], r_perc)
}

# same as plot_map_20panel but first convert difference
# layers into percent change
plot_map_20panel_perc <- function(
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
    midpoint_diff = 0
) {
  
  r_perc <- calc_rast_perc_change(r = r, info = info, 
                                  type_absolute = type_absolute, 
                                  type_diff = type_diff)
  
  if(is.null(name4diff)) {
    name4diff <- paste("\u0394", name4absolute, '(%)')
  }
  
  plot_map_20panel_perc(
    r = r_perc,
    info = info,
    type_absolute = type_absolute,
    type_diff = type_diff,
    title = title,
    name4absolute = name4absolute,
    name4diff = name4diff,
    legend_title_absolute = legend_title_absolute,
    legend_title_diff = legend_title_diff,
    palette_absolute = palette_absolute,
    palette_diff = palette_diff,
    lims_diff = lims_diff,
    midpoint_diff = midpoint_diff
  )
}


# four panels, top row is is historical c3, bottomr row is future c9,
# columns are two grazing levels
plot_4panel_c3c9 <- function(r_c3, r_c9, info_c3, info_c9) {
  
  stopifnot(
    nrow(info_c3) == 2,
    nrow(info_c9) ==2
  )
  
  info_c3$type <- 'c3'
  info_c9$type <- 'c9'
  info_comb <- bind_rows(info_c3, info_c9) %>% 
    arrange(type, RCP, years, graze) %>%
    ungroup() %>% 
    mutate(tag_label = fig_letters[1:n()])
  
  levs <- levels(info_comb$graze)
  graze_labs <- paste0(levs[levs %in% info_comb$graze], '\nGrazing')
  names4leftside <- c('c3' = 'SEI class',
                      'c9' = 'Change in SEI class')
  levs_l <- list(c3 = data.frame(ID = 1:3, 
                                 class = c3Names),
                 c9 = data.frame(ID = 1:9, 
                                 class = c9Names))
  label_left <- info_comb %>% 
    filter(.data$graze == unique(.data$graze)[1]) %>% 
    mutate(rcp_year = rcp_label(.data$RCP, .data$years,
                                include_parenth = FALSE),
           label_left = paste0(rcp_year, '\n', names4leftside[type])) %>% 
    pull(label_left)
  
  fill <- list(
    'c3' = scale_fill_c3,
    'c9' = scale_fill_c9
  )
  
  r_l <- list('c3' = r_c3,
              'c9' = r_c9)
  
  guide_l <- list('c3' = 'right', 'c9' = 'none')
  
  plots1 <- pmap(info_comb[c('id', 'type', 'tag_label')], 
                 function(id, type, tag_label) {
                   fl <- fill[[type]]
                   r <- r_l[[type]][[id]]
                   levels(r) <- levs_l[[type]]
                   plot_map2(r = r,
                             panel_tag = tag_label) + 
                     fl() + 
                     theme(legend.position = guide_l[[type]],
                           plot.tag.position =  'topleft',
                           plot.tag.location = 'panel',
                           plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
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
  color_matrix1 <- color_matrix() # to serve as a legend
  # this is hard coded, would need to make this more
  # flexible if want to change number of panels
  s <- list(plot_spacer())
  
  plots2 <- c(s, plots_top, s,
              plots_left[1], plots1[1:2], list(guide_area()), 
              plots_left[2], plots1[3:4], list(free(color_matrix1)))
  
  # for troubleshooting layout
  # plots2 <- map(plots2, \(x) x +   theme(plot.background = element_rect(color = "blue", fill = NA, size = 1)))

  patchwork::wrap_plots(plots2, ncol = 4,
                             widths = c(0.13, 1, 1, 0.55), 
                             heights = c(0.13, 1, 1)) + 
    plot_layout(guides = 'collect')
 
}



# within a climate scenario comparisons to a grazing reference level (gref)
# top row is is historical, bottom is future
# left column is c3 moderate grazing, other columns (can be 1 or more,
# depending on what 'target' grazing grazing levels are selected
plot_c3c9_gref <- function(r_c3, r_c9, info_c3, info_c9) {
  
  stopifnot(
    nrow(info_c3) == 2
  )
  
  info_c3$type <- 'c3'
  info_c9$type <- 'c9'
  info_comb <- bind_rows(info_c3, info_c9) %>% 
    arrange(RCP, years, type, graze) %>%
    ungroup() %>% 
    mutate(tag_label = fig_letters[1:n()])
  
  graze_labs <- paste0(unique(info_comb$graze), ' Grazing')
  n_target <- lu(info_c9$graze)
  names4top <- paste0(graze_labs, 
                      paste0('\n',c('(SEI class)', 
                                    rep('(Change in SEI class)', n_target)))
  )
  levs_l <- list(c3 = data.frame(ID = 1:3, 
                                 class = c3Names),
                 c9 = data.frame(ID = 1:9, 
                                 class = c9Names))
  label_left <- info_comb %>% 
    filter(.data$graze == unique(.data$graze)[1]) %>% 
    mutate(rcp_year = rcp_label(.data$RCP, .data$years,
                                include_parenth = FALSE)) %>%  
    pull(rcp_year)
  
  fill <- list(
    'c3' = scale_fill_c3,
    'c9' = scale_fill_c9
  )
  
  r_l <- list('c3' = r_c3,
              'c9' = r_c9)
  
  guide_l <- list('c3' = 'right', 'c9' = 'none')
  
  plots1 <- pmap(info_comb[c('id', 'type', 'tag_label')], 
                 function(id, type, tag_label) {
                   fl <- fill[[type]]
                   r <- r_l[[type]][[id]]
                   levels(r) <- levs_l[[type]]
                   plot_map2(r = r,
                             panel_tag = tag_label) + 
                     fl() + 
                     theme(legend.position = guide_l[[type]],
                           plot.tag.position =  'topleft',
                           plot.tag.location = 'panel',
                           plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
                 })
  
  # labels for the top margin
  plots_top <- map(names4top, function(x) {
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
  color_matrix1 <- color_matrix() # to serve as a legend
  # this is hard coded, would need to make this more
  # flexible if want to change number of panels
  s <- list(plot_spacer())
  
  n <- length(plots1)
  if(n == 4) {
    plots2 <- c(s, plots_top,
                plots_left[1], plots1[1:(n/2)], 
                plots_left[2], plots1[(n/2 + 1):n], 
                s, list(guide_area()), list(free(color_matrix1)))
  } else {
    plots2 <- c(s, plots_top,
                plots_left[1], plots1[1:(n/2)], 
                plots_left[2], plots1[(n/2 + 1):n], 
                s, list(guide_area()), list(free(color_matrix1)), rep(s, n/2 - 2))
  }
  
  patchwork::wrap_plots(plots2, ncol = 1 + n/2,
                        widths = c(0.13, rep(1, n/2)), 
                        heights = c(0.13, 1, 1, 0.6)) + 
    plot_layout(guides = 'collect')
  
}

# crs -----------------------------------------------------------------------

# the crs to be used for the sagebrush conservation design (this is the same
# one as used by NLCD)
crs_scd <- terra::crs("PROJCRS[\"Albers_Conical_Equal_Area\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Albers Equal Area\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",23,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-96,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")



# polygons for basemaps etc---------------------------------------------------

states <- sf::st_transform(sf::st_as_sf(spData::us_states), crs = crs_scd)

load_template <- function(template_path = file.path(
  "data_processed",  "interpolated_rasters","biomass", "v4",
  "fire1_eind1_c4grass1_co20_2503v4_Sagebrush_biomass_Current_Current_Light_Current.tif")) {
  # any of the correct version of biomass interpolated layers works
  # as a template
  r <- terra::rast(template_path)[[1]]
  r[!is.na(r)] <- 1
  r
}

load_wafwa_ecoregions <- function(total_region = FALSE, wafwa_only = FALSE) {
  # file created by the 01_ecoregions.R script
  shp1 <- sf::st_read('data_processed/ecoregions/four_regions_v1.gpkg')
  shp2 <- shp1|> 
    dplyr::rename(geometry = geom)
  
  if(wafwa_only) {
    shp2$ecoregion <- shp2$wafwa_NAME
    shp2 <- shp2 |> 
      dplyr::group_by(ecoregion) |> 
      dplyr::summarize(geometry = sf::st_union(.data$geometry),
                       .groups = 'drop')
  } else {
    shp2$ecoregion <- shp2$region
  }


    # also have a seperate row for the entire region
  if(total_region) {
    shp2 <- shp2 |> 
      dplyr::summarize(geometry = sf::st_union(.data$geometry)) |> 
      dplyr::mutate(ecoregion = "Entire study area") |> 
      dplyr::bind_rows(shp2) |> 
      dplyr::mutate(ecoregion = region_factor(as.character(.data$ecoregion)))
  } else {
    # these factor levels not being alphabetical could cause problems
    # in some code
    shp2$ecoregion <- region_factor(shp2$ecoregion, include_entire = FALSE,
                                    wafwa_only = wafwa_only) 
  }
  shp2
}

load_wafwa_ecoregions_raster <- function(
    wafwa_only = FALSE) {
  
  r_template <- load_template()
  eco1 <- load_wafwa_ecoregions(wafwa_only = wafwa_only)
  r_eco <- terra::rasterize(
    terra::vect(eco1),
    r_template,
    field = 'ecoregion',
    touches = FALSE
  )
  
  r_eco[is.na(r_template)] <- NA
  
  # reordering so that levels
  # in the raster are the same order as the factor
  # from region_factor
  lvls_old <- as.data.frame(levels(r_eco)[[1]])
  desired_levs <- levels(eco1$ecoregion)
  order <- 0:(length(desired_levs) -1)
  names(order) <- desired_levs
 
  lvls_old$ID_to <- order[lvls_old$ecoregion]
  lvls_old$ecoregion_to <- desired_levs[lvls_old$ID_to + 1]

  r_eco2 <- terra::classify(r_eco,
                            as.matrix(lvls_old[c('ID', 'ID_to')]))
  
  lvls <- lvls_old[c('ID_to', 'ecoregion_to')]
  names(lvls) <- stringr::str_replace(names(lvls), "_to", "")
  lvls <- dplyr::arrange(lvls, .data$ID)

  levels(r_eco2) <- lvls
  r_eco2
}
load_cell_size <- function() {
  r <- load_wafwa_ecoregions_raster()
  s <- terra::cellSize(r, mask = TRUE,
                       unit = 'ha',
                       transform = FALSE)
  s
}
