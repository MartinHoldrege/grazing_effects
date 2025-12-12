# functions for extracting summary data from rasters
# and for plotting that data to create figures of the
# distributions (e.g. histograms)

# get data needed for plotting a histogram
hist_from_rast <- function(r, breaks, by_ecoregion = TRUE) {
  stopifnot(inherits(r, "SpatRaster"))
  
  # computation is 1 layer at a time for memory reasons
  nms <- names(r)
  
  # Identify ecoregion layer & value layers
  if (by_ecoregion) {
    stopifnot("ecoregion" %in% nms)
    eco_layer <- r[["ecoregion"]]
    val_names <- setdiff(nms, "ecoregion")
  } else {
    val_names <- nms
  }
  
  # Read ecoregion values once
  if (by_ecoregion) {
    eco_vec <- terra::values(eco_layer, mat = FALSE)
    # drop NAs; assume numeric / factor-like IDs
    regions <- sort(unique(eco_vec))
    
    levs <- levels(eco_layer)[[1]]
    regions <- levs$ID
    
    # Set up output list: one element per region + "Entire study area"
    region_labels <- c(levs$ecoregion, "Entire study area")
    out <- setNames(vector("list", length(region_labels)), region_labels)
    for (nm_reg in region_labels) {
      out[[nm_reg]] <- setNames(vector("list", length(val_names)), val_names)
    }
    
    # Loop over value layers 
    for (nm in val_names) {
      v <- terra::values(r[[nm]], mat = FALSE)
      
      # Entire study area
      out[["Entire study area"]][[nm]] <- hist(v, breaks = breaks, plot = FALSE)
      
      # By region
      for (reg in regions) {
        idx <- eco_vec == reg
        # guard against regions with all NA in this layer
        x_reg <- v[idx]
        reg_name <- levs$ecoregion[levs$ID == reg]
        if (all(is.na(x_reg))) {
          out[[reg_name]][[nm]] <- hist(numeric(0), breaks = breaks, plot = FALSE)
        } else {
          out[[reg_name]][[nm]] <- hist(x_reg, breaks = breaks, plot = FALSE)
        }
      }
    }
    
  } else {
    # No ecoregion grouping: just loop over layers 
    out <- setNames(vector("list", length(val_names)), val_names)
    for (nm in val_names) {
      v <- terra::values(r[[nm]], mat = FALSE)
      out[[nm]] <- hist(v, breaks = breaks, plot = FALSE)
    }
  }
  
  out
}


# convert nested list output from hist_from_rast
# to a dataframe
hist_dat2df <- function(l, by_ecoregion = TRUE,
                        into = c('type', 'RCP', 'years', 'graze', 'GCM')) {
  
  inner2df <- function(x) {
    stopifnot(c('counts', 'density', 'mids') %in% names(x))
    tibble(counts = x$counts,
           density = x$density,
           mids = x$mids)
  }
  
  id2df <- function(l) {
    map_dfr(l, .f = inner2df, .id = 'id')
  }
  
  region2df <- function(l) {
    map_dfr(l, .f = id2df, .id = 'region')
  }
  
  f <- if(by_ecoregion) region2df else id2df
  df <- f(l)
  
  out <- create_rast_info(unique(out$id), 
                   into = into) %>% 
    right_join(df, by = 'id')
  out
}

