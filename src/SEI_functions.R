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
            nlyr(ecor_raster) == 1)
  
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


# continuous sei to three categories of SEI
sei2c3 <- function(x) {
  stopifnot(
    x>=0,
    x<=1
  )
  
  c3 <- dplyr::case_when(
    x > 0.431 ~ 'CSA',
    x <= 0.431 & x > 0.173 ~ 'GOA',
    x <= 0.173 ~ 'ORA'
  )
  
  factor(c3, levels = c('CSA', "GOA", "ORA"))
  
}





