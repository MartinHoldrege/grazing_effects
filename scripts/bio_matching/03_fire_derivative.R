# Purpose: See what effect quantile matching has (spatially) on 
# the derivative of the fire equation for both perennials and annuals


# Author: Martin Holdrege

# Started: November 11, 2024


# parameters --------------------------------------------------------------

run <- "fire1_eind1_c4grass1_co20_2311"
version_interp <- 'v3'
graze <- c('grazL' = 'Light')
run2 <- paste0(run, "_", names(graze))

mask_cutoff <- 1 # matching quality cutoff used for the mask
# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(patchwork)
source('src/fire_functions.R')
source("src/mapping_functions.R")
source("src/fig_params.R")

# read in data ------------------------------------------------------------

# functions for quantile matching
path <- 'data_processed/temp_rds/qm_funs.rds'
if(file.exists(path)) {
  qm_l <- readRDS(path)
} else {
  source('scripts/bio_matching/02_quantile_matching.R')
  qm_l <- readRDS(path)
}

qual_cutoff <- qm_l$qual_cutoff
rcp <- 'Current'

# interpolated climate 
vars_clim <- c("MAT", "MAP", "psp")
paths_clim <- file.path('data_processed/interpolated_rasters/climate', 
                   paste0(vars_clim, 
                          "_climate_", rcp, "_Current_Current_20240605_",
                          version_interp, ".tif"))
clim1 <- rast(paths_clim)
names(clim1) <- tolower(vars_clim)

# interpolated biomass
pfts <- c("pfg" = "Pherb", "afg" = "Aherb")
paths_pft <- file.path("data_processed/interpolated_rasters/biomass", 
                        paste0(run, version_interp, '_', pfts, 
                               "_biomass_", rcp, "_Current_",
                               graze, "_Current.tif"))
bio1 <- rast(paths_pft)
names(bio1) <- names(pfts) # pfg/afg naming needed for functions used below

pred_raw <- c(clim1, bio1) # fire equation predictor variables

# for the maskmask
r_qual <- rast(paste0("data_processed/interpolation_quality/matching_quality",
                      version_interp, ".tif"))


# functions ---------------------------------------------------------------

mm <- function(r) as.numeric(minmax(r))

# create mask -------------------------------------------------------------

mask <- r_qual
mask[mask > mask_cutoff] <- NA
mask[mask <= mask_cutoff] <- 1

# quantile mapping --------------------------------------------------------
qm_Pherb <- qm_l$Pherb
qm_Aherb <- qm_l$Aherb

pred_raw <- trim(mask(pred_raw, mask))
pred_qm <- pred_raw # will contain the quantile mapped biomass values
pred_qm[['pfg']][] <- qm_Pherb(values(pred_raw$pfg))
pred_qm[['afg']][] <- qm_Aherb(values(pred_raw$afg))


# fire probability --------------------------------------------------------

# 'raw' predictor variables
fire_prob_raw <- fire_rast_wrapper(f = predict_fire,
                                   r = pred_raw,
                                   run_checks = FALSE)*100

fire_prob_qm <- fire_rast_wrapper(f = predict_fire,
                                  r = pred_qm,
                                  run_checks = FALSE)*100


# derivative of fire probability------------------------------------------
pred_vars <- names(pred_raw)
names(pred_vars) <- pred_vars
# for naming in figures
pred_var_names <- c(toupper(pred_vars[1:3]), pfts[pred_vars[4:5]])
names(pred_var_names) <- pred_vars

# percentage point change in fire probability with 1 unit change in deriv_var
# (except for psp, it's per 0.1 unit)
fire_deriv_raw <- map(pred_vars, function(deriv_var) {
  p <- fire_rast_wrapper(f = calc_fire_deriv,
                    r = pred_raw,
                    deriv_var = deriv_var,
                    run_checks = FALSE)
  p*100 
}) %>% 
  rast()

# derivatives calculated using quantile mapped biomass
fire_deriv_qm <- map(pred_vars, function(deriv_var) {
  p <- fire_rast_wrapper(f = calc_fire_deriv,
                         r = pred_qm,
                         deriv_var = deriv_var,
                         run_checks = FALSE)
  p*100 
}) %>% 
  rast()

# figure components -------------------------------------------------------

fill_slope <- function(limits = NULL, oob = scales::censor) {
  scale_fill_gradient2(low = "blue",
                       high = 'red',
                       mid = 'lightgrey',
                       midpoint = 0,
                       na.value = 'transparent',
                       name = 'Slope',
                       limits = limits,
                       oob = oob) 
}

smaller_legend <- function() {
  theme(
    legend.title = element_text(size = rel(0.5)),  # Reduce legend title size
    legend.text = element_text(size = rel(0.5)),   # Reduce legend text size
    legend.key.size = unit(0.5, "lines"),          # Reduce legend key size (relative to lines)
    legend.spacing = unit(0.5, "lines")           # Reduce spacing between legend items
  )
}

base_deriv <- function(limits = NULL, oob = scales::censor) {
  list(
    fill_slope(limits = limits, oob = oob),
    smaller_legend(),
    theme(legend.position = 'right')
  )
}



# maps predictor vars -----------------------------------------------------
lims_pred <- map(pred_vars, function(var) {
  range(c(mm(pred_raw[[var]]), 
          mm(pred_qm[[var]]))
  )
})

base_pred <- function(var) {
  list(smaller_legend(),
       scale_fill_gradient(low = cols_pred[[var]]['low'],
                        high = cols_pred[[var]]['high'],
                        na.value = 'transparent',
                        name = pred_var_names[var],
                        limits = lims_pred[[var]]),
       theme(legend.position = 'right'))
}

maps_pred_raw <- map(pred_vars, function(var) {
  plot_map2(pred_raw[[var]],
            maintitle = pred_var_names[var]) +
    base_pred(var)

})

maps_pred_qm <- map(pred_vars, function(var) {
  plot_map2(pred_qm[[var]],
            maintitle = pred_var_names[var]) +
    base_pred(var)
  
})

cap0 <- paste('data from', rcp, 
      'conditions and', run2,
      '\nStudy area = SCD cells with <', mask_cutoff, 'matching quality')

pdf(paste0('figures/bio_matching/fire_pred_vars_qm', qual_cutoff, '.pdf'),
    width = 10, height = 7)
g1 <- wrap_plots(maps_pred_raw, nrow = 2)
g1 + plot_annotation(title = 'Fire predictor variables',
                     subtitle = '(Raw STEPWAT2 biomass)',
                     caption = paste(cap0, '\n\n'))

g2 <- wrap_plots(maps_pred_qm, nrow = 2)
g2 + plot_annotation(title = 'Fire predictor variables',
                     subtitle = '(Quantile mapped STEPWAT2 biomass)',
                     caption = paste(cap0, '\n', qm_l$caption))
dev.off()

# maps--derivatives -------------------------------------------------------

lims_deriv <- map(pred_vars, function(var) {
  range(c(mm(fire_deriv_raw[[var]]), 
             mm(fire_deriv_qm[[var]]))
        )
})

maps_deriv_raw <- map(pred_vars, function(var) {
  g1 <- plot_map2(fire_deriv_raw[[var]],
            maintitle = pred_var_names[var]) +
    # setting constant limits so maps are comparable
    # between qm and raw datasets for a given predictor var
    base_deriv(limits = lims_deriv[[var]])
  
  g2 <- plot_map2(fire_deriv_raw[[var]],
                  maintitle = pred_var_names[var]) +
    base_deriv(limits = c(-0.1, 0.1), oob = scales::squish)
  
  # variable (among predictors) and fixed limits
  list(variable = g1, fixed = g2)
})

maps_deriv_qm <- map(pred_vars, function(var) {
  g1 <- plot_map2(fire_deriv_qm[[var]],
                  maintitle = pred_var_names[var]) +
    # setting constant limits so maps are comparable
    # between qm and raw datasets for a given predictor var
    base_deriv(limits = lims_deriv[[var]])
  
  g2 <- plot_map2(fire_deriv_qm[[var]],
                  maintitle = pred_var_names[var]) +
    base_deriv(limits = c(-0.1, 0.1), oob = scales::squish)
  
  # variable (among predictors) and fixed limits
  list(variable = g1, fixed = g2)
})

cap1 <- paste('predictions made using simulation data from', rcp, 
              'conditions', run2,
              '\nStudy area = SCD cells with <', mask_cutoff, 'matching quality')

annotation_raw <- function(prefix = '') {
  plot_annotation(
  title = 'Derivative of fire equation (raw stepwat biomass used)',
  subtitle = 'change in fire probability (delta # fires/100 years) per 1% change in predictor variable',
  caption = paste(prefix, cap1, '\n\n'))
}

annotation_qm <- function(prefix = '') {
  plot_annotation(
  title = 'Derivative of fire equation (quantile mapped stepwat biomass used)',
  subtitle = 'change in fire probability (delta # fires/100 years) per 1% change in predictor variable',
  caption = paste(prefix, cap1, 
                  '\n', qm_l$caption)) 
}

pdf(paste0('figures/bio_matching/fire_slopes_qm', qual_cutoff, '.pdf'),
    width = 10, height = 7)
g1 <- wrap_plots(map(maps_deriv_raw, \(x) x$fixed), nrow = 2)
g1 + annotation_raw(prefix = 'Fixed color scales')

g2 <- wrap_plots(map(maps_deriv_qm, \(x) x$fixed), nrow = 2)
g2 + annotation_qm(prefix = 'Fixed color scales')

g1 <- wrap_plots(map(maps_deriv_raw, \(x) x$variable), nrow = 2)
g1 + annotation_raw()

g2 <- wrap_plots(map(maps_deriv_qm, \(x) x$variable), nrow = 2)
g2 + annotation_qm()

dev.off()

