# Purpose: fit ecdfs to RAP and stepwat biomass so can do quantile matching

# Martin Holdrege

# Started: Nov 8, 2024

# params ------------------------------------------------------------------

dataset <- 'interpolated'

# dependencies ------------------------------------------------------------

library(tidyverse)

# read in data ------------------------------------------------------------

# list created in 01_combine_data.R
m <- readRDS('data_processed/temp_rds/rap_sw_matching.rds') 

# vectors -----------------------------------------------------------------

pfts <- c('Pherb' = 'Pherb', 'Aherb' = 'Aherb')

# fitting cdfs ------------------------------------------------------------

# rap datasets
rap_df <- m$rap_comb %>% 
  filter(dataset == !!dataset)

sw_df <- m$sw_comb %>% 
  filter(dataset == !!dataset)

# fit emperical cdfs
cdfs <- list()
cdfs$rap <- map(pfts, function(pft) {
  x <- rap_df[rap_df$PFT == pft, ]$biomass
  ecdf(x)
})

cdfs$sw <- map(pfts, function(pft) {
  x <- sw_df[sw_df$PFT == pft, ]$biomass
  ecdf(x)
})

# functions ---------------------------------------------------------------

qm_factory <- function(pft) {
  # using this approach so the ecdfs are encapsulated in the function environment
  cdf_from <- cdfs$sw[[pft]] 
  cdf_to <- cdfs$rap[[pft]]
  function(x) {
    probs <- cdf_from(x)
    values <- quantile(cdf_to, probs = probs)
    as.numeric(values)
  }
}

# create functions where input is a vector of biomass
# and output is a vector of adjusted biomass, so that it
# is on the same scale as RAP
qm_Aherb <- qm_factory('Aherb')
qm_Pherb <- qm_factory('Pherb')

# save objects ------------------------------------------------------------

out <- list('Aherb' = qm_Aherb,
            'Pherb' = qm_Pherb,
            'dataset' = dataset)
saveRDS(out, 'data_processed/temp_rds/qm_funs.rds')
