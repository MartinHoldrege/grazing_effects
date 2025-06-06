source('src/qm_functions.R')


# the following list was created by 
# running the parse_q_curves() function in the SEI/src/general_functions.R
# script, and then running dput(nameOfList) (and change Aherb and Pherb names)
q_curves_cover <- list(Sagebrush = structure(list(cover = c(0, 0.02, 0.03, 0.05, 0.1, 
0.15, 0.2, 0.25, 0.5, 0.75, 1), great_basin = c(0, 0.05, 0.1, 
0.5, 1, 1, 1, 1, 1, 1, 1), intermountain = c(0, 0.05, 0.1, 0.5, 
0.75, 0.9, 0.95, 1, 1, 1, 1), great_plains = c(0.33, 0.5, 0.7, 
0.9, 0.95, 1, 1, 1, 1, 1, 1)), row.names = c(NA, -11L), class = "data.frame"), 
    Pherb = structure(list(cover = c(0, 0.02, 0.03, 0.05, 0.1, 
    0.15, 0.2, 0.25, 0.5, 0.75, 1), great_basin = c(0, 0.02, 
    0.05, 0.45, 1, 1, 1, 1, 1, 1, 1), intermountain = c(0, 0.02, 
    0.05, 0.45, 0.7, 0.88, 0.94, 1, 1, 1, 1), great_plains = c(0, 
    0, 0, 0.2, 0.7, 0.88, 0.94, 0.95, 1, 1, 1)), row.names = c(NA, 
    -11L), class = "data.frame"), Aherb = structure(list(cover = c(0, 
    0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 1), great_basin = c(1, 
    0.98, 0.95, 0.9, 0.5, 0.25, 0.1, 0, 0, 0, 0), intermountain = c(1, 
    0.98, 0.95, 0.9, 0.5, 0.25, 0.1, 0, 0, 0, 0), great_plains = c(1, 
    1, 1, 1, 1, 0.75, 0.55, 0.4, 0.2, 0, 0)), row.names = c(NA, 
    -11L), class = "data.frame"))

q_curves_cover2 <- map(q_curves_cover, function(x) {
  x %>% rename("Great Plains" = 'great_plains',
               "Intermountain West" = 'intermountain', 
                "Southern Great Basin" = 'great_basin')
})

bio2cov_Aherb <- function(x) {
  bio_qm <- qm_Aherb(x)
  # from equation in Holdrege et al 2024 REM (appendix C)
  cover <- 2.457+0.6319*bio_qm
  cover
}

bio2cov_Pherb <- function(x) {
  bio_qm <- qm_Pherb(x)
  # from equation in Holdrege et al 2024 REM (appendix C)
  cover <- 4.833 + 0.5073*bio_qm
  cover
}

# biomass to cover from Holdrege et al 2024 REM (appendix C)
bio2cov_hold <- function(x, pft) {
  betas <- list(
    'Aherb' = c('B0' = 2.457, 'B1' = 0.6319),
    'Pherb' = c('B0' = 4.833, 'B1' = 0.5073),
    'Sagebrush' = c('B0' = 0.1078, 'B1' = 0.03451)
  )
  stopifnot(pft %in% names(betas))
  b <- betas[[pft]]
  b['B0'] + b['B1']*x
}


# convert biomass to cover, including quantile mapping of the
# biomass to the RAP/RCMAP scale (not for sagebrush
# it's a straight quantile mapping to cover, with the others
# the quantile mapped biomass gets converted to cover
# via a linear equation)
bio2cov <- function(x, pft) {
  
  stopifnot(
    pft %in% c("Aherb", "Pherb", "Sagebrush")
  )
  dplyr::case_when(
    pft == 'Aherb' ~ bio2cov_Aherb(x),
    pft == 'Pherb' ~ bio2cov_Pherb(x),
    pft == 'Sagebrush' ~ qm_Sagebrush_bio2cov(x),
  )
}

#' convert stepwat biomass to Q (quality) score
#'
#' @param x biomass
#' @param region ecoregion
#' @param pft string
bio2q <- function(x, pft, region) {
  
  tab <- q_curves_cover2[[pft]]
  
  region <- region2wafwa(region)
  
  cover <- bio2cov(x, pft)/100 # convert from % to proportion
  stopifnot(region %in% names(tab),
            is.character(pft),
            is.character(region),
            length(pft) == 1,
            length(region) == 1)
  tab$cover
  # calculate q values via linear interpolation
  q <- qm_from_quantiles(cover, 
                         from = tab$cover, # cover values
                         to = tab[[region]],  # q values in table
                         ascending = FALSE)
  
  q
}

#' convert cover to Q (quality) score
#'
#' @param x cover
#' @param region ecoregion
#' @param pft string
cov2q <- function(x, pft, region) {
  
  tab <- q_curves_cover2[[pft]]
  
  region <- region2wafwa(region)
  
  stopifnot(region %in% names(tab),
            is.character(pft),
            is.character(region),
            length(pft) == 1,
            length(region) == 1)
  # calculate q values via linear interpolation
  q <- qm_from_quantiles_fast(x, 
                         from = tab$cover, # cover values
                         to = tab[[region]],  # q values in table
                         ascending = FALSE)
  
  q
}

#' calculate Q and SEI from rasters of biomass
#'
#' @param r raster with 3 layers (biomass of sage, annuals and perennials)
#' @param eco_raster raster providing the 3 ecoregions
#' @param type_string part of layer name giving the data category
#'
#' @returns
#' 4 layer raster (3 Qs and SEI)
bio2qsei_raster <- function(r, eco_raster, type_string = 'biomass') {
  stopifnot(isTRUE(terra::compareGeom(r, eco_raster)))
  
  pfts <- c('Sagebrush', 'Pherb', 'Aherb')
  
  regex <- paste(pfts, collapse = '|')
  
  nms <- names(r)
  stopifnot(str_detect(nms, regex),
            str_detect(nms, type_string),
            nlyr(r) == 3,
            nlyr(eco_raster) == 1)
  
  regions <- levels(eco_raster)[[1]]$ecoregion
  stopifnot(length(regions) == 3)
  
  comb <- tidyr::expand_grid(
    PFT = pfts,
    region = regions
    
  )
  
  r_q <-  r
  r_q[] <- NA
  
  for (i in 1:nrow(comb)) {
    pft <- comb$PFT[i]
    region <- comb$region[i]
    lyr <- str_subset(nms, pattern = pft)
    
    biomass <- as.vector(r[[lyr]][eco_raster == region])
    
    r_q[[lyr]][eco_raster == region] <- bio2q(biomass, pft = pft, 
                                              region = region)
    
  }
  
  stopifnot(
    nlyr(r_q) == 3,
    names(r_q) == nms
  )
  r_sei <- app(r_q, fun = prod) # multiply 3 q's together
  
  # renaming file name components
  names(r_q) <- str_replace(nms, type_string, 'Q')
  
  # the PFT part of the layer name becomes NA for SEI
  names(r_sei) <- str_replace(nms[1], regex, 'NA') %>% 
    str_replace(type_string, "SEI")
  
  c(r_q, r_sei)
}


#' convert cover to quality score
#'
#' @param r raster of cover (can be multiple bands, but all must be for one pft) 
#' @param pft the plant functional group that the cover is for
#' @param eco_raster raster providing ecoregions
cov2q_raster <- function(r, eco_raster, pft) {
  stopifnot(isTRUE(terra::compareGeom(r, eco_raster)))
  
  regions <- levels(eco_raster)[[1]]$ecoregion
  stopifnot(length(regions) == 3)
  
  ids <- names(r)
  result_matrix <- matrix(NA, nrow = ncell(r), ncol = nlyr(r))
  colnames(result_matrix) <- ids
  
  for (region in regions) {
    mask <- as.vector(eco_raster == region)
    mask[is.na(mask)] <- FALSE
    tmp <- r
    tmp[!mask] <- NA
    result <- app(tmp, cov2q, pft = pft, region = region)
    result_matrix[mask, ] <- as.matrix(result)[mask, ]
  }
  
  # Fill raster all at once
  r_q <- r
  r_q[] <- NA
  values(r_q) <- result_matrix
  r_q
}


# continuous sei to three categories of SEI
sei2c3 <- function(x) {
  
  if(is.numeric(x)) {
    stopifnot(
      (x>=0 & x<=1) | is.na(x)
    )
    
    c3 <- dplyr::case_when(
      x > 0.431 ~ 'CSA',
      x <= 0.431 & x > 0.173 ~ 'GOA',
      x <= 0.173 ~ 'ORA',
      TRUE ~ NA
    )
    
    c3 <- factor(c3, levels = c('CSA', "GOA", "ORA"))
    
  } else if("SpatRaster" %in% class(x)) {
    m <- matrix(c(
      -Inf, 0.173, 3,  # ORA
      0.173, 0.431, 2, # GOA
      0.431, Inf, 1    # CSA
    ), ncol = 3, byrow = TRUE)
    
    c3 <- terra::classify(x, m, right = TRUE)
    
    # levels(c3) <- data.frame(
    #   ID = 1:3,
    #   c3 = c("CSA", "GOA", "ORA")
    # )
    
    names(c3) <- names(x)

  } else {
    stop('class not supported')
  }
  
  c3
}

c3_factor <- function(x) {
  levels <- c('CSA', "GOA", "ORA")
  out <- if(all(x %in% c(1:3))) {
    factor(x, levels = 1:3, labels = levels)
  } else if (all(x %in% levels)) {
    factor(x, levels)
  } else {
    stop('incorrect c3 values')
  }
  out
}



# takes two SEI class rasters, and creates a 9 class change raster
c3toc9 <- function(current, future) {
  c9From <-  c(11, 12, 13, 21, 22, 23, 31, 32, 33); 
  c9To <-  c(1, 2, 3, 4, 5, 6, 7, 8, 9);
  
  current10 <- current*10
  
  c9a <- current10 + future
  
  rcl <- matrix(c(c9From, c9To), byrow = FALSE, ncol = 2)
  c9b <- terra::classify(c9a, rcl)
  # c9Names <-  c(
  #   'Stable CSA',
  #   'CSA becomes GOA',
  #   'CSA becomes ORA',
  #   'GOA becomes CSA',
  #   'Stable GOA',
  #   'GOA becomes ORA',
  #   'ORA becomes CSA',
  #   'ORA becomes GOA',
  #   'Stable ORA'
  # )
  # 
  # levels(c9b) <- data.frame(
  #   ID = 1:3,
  #   c9 = c9Names
  # )
  names(c9b) <- names(future)
  c9b
}

# convert two digits into c3 (first digit is c3)
c3eco_to_c3 <- function(x) {
  stopifnot(x >= 10 & x < 40)
  c3_num <- round(x/10, 0) 
  c3_factor(c3_num)
}

# convert two digit number into the region
# need to make sure the regions vector is from
# the raster used for the zonal stats, so the order is correct
c3eco_to_eco <- function(x, regions) {
  stopifnot(x >= 10 & x < 40)
  eco_num <- x %% 10 # returning the digit in the ones position
  stopifnot(sort(unique(eco_num)) == 0:(length(regions) - 1))
  # the IDs in factor spatrasters start at 0, not 1
  regions[eco_num + 1]
}

c9_factor <- function(x) {
  c9Names <-  c(
    'Stable CSA',
    'CSA becomes GOA',
    'CSA becomes ORA',
    'GOA becomes CSA',
    'Stable GOA',
    'GOA becomes ORA',
    'ORA becomes CSA',
    'ORA becomes GOA',
    'Stable ORA'
  )
  
  if(all(x %in% c9Names)) {
    levels = c9Names
  } else if (all(x %in% 1:9)) {
    levels <- 1:9
  } else {
    stop('x does not have the right values')
  }
  
  factor(x, levels = levels,
         labels = c9Names)
}

percent_csa <- function(x, na.rm = TRUE) {
  mean(x > 0.431, na.rm = na.rm)*100
}

percent_goa <- function(x, na.rm = TRUE) {
  mean(x > 0.173 & x <= 0.431, na.rm = na.rm)*100
}


load_scd_cover <- function() {
  # file created in 04_scd_aggregate.R
  terra::rast('data_processed/scd/SEIv11_2017_2020_1000_20211228_cover560.tif')
}

load_scd_q <- function() {
  # file created in 04_scd_aggregate.R
  # bands:
  # SEI
  # Q1 --sage
  # Q2 --perennials
  # Q3 --annuals
  # Q4 --human mod
  # Q5 --trees
  terra::rast('data_processed/scd/SEI-Q_v11_2017-2020_1000.tif')
}


