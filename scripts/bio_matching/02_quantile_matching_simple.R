# Purpose: create output for use in STEPWAT for for quantile mapping
# based on using just a select few quantiles. 
# see the 02_quantile_matching.R script for more extensive 
# exploration

# Martin Holdrege

# Started: Dec 3, 2024

# params ------------------------------------------------------------------

# probabilities to calculate the quantiles of
probs <- c(0, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98, 0.995, 1)
qual_cutoff <- 0.5 # matching quality cutoff used determine how good
# the interpolation needed to be

# dependencies ------------------------------------------------------------

library(tidyverse)
source('src/qm_functions.R')

# read in data ------------------------------------------------------------

# files created in 01_combine_data.R
sw1 <- read_csv(paste0('data_processed/qm/stepwat_for_qm_', qual_cutoff, 
                     'match.csv'),
                show_col_types = FALSE)

rap1 <- read_csv(paste0('data_processed/qm/rap_for_qm_', qual_cutoff, 
                      'match.csv'),
                 show_col_types = FALSE)

# vectors -----------------------------------------------------------------

pfts <- c('Pherb' = 'Pherb', 'Aherb' = 'Aherb')

# fitting cdfs ------------------------------------------------------------

# fit emperical cdfs to interpolated data
cdfs <- list()

cdfs$rap <- map(pfts, function(pft) {
  x <- rap1[rap1$PFT == pft, ]$biomass
  ecdf(x)
})

cdfs$sw <- map(pfts, function(pft) {
  x <- sw1[sw1$PFT == pft, ]$biomass
  ecdf(x)
})

# quantiles
quants <- map(cdfs, function(l) {
  map(l, \(x) as.numeric(quantile(x, probs = probs)))
})

# convert to a dataframe for outputting
to_long <- function(l, values_to) {
  l %>% 
    bind_cols() %>% 
    mutate(quantile = probs) %>% 
    pivot_longer(cols = all_of(unname(pfts)),
                 names_to = 'PFT',
                 values_to = values_to)
}

out <- to_long(quants$rap, values_to = 'biomass_rap') %>% 
  left_join(to_long(quants$sw, values_to = 'biomass_stepwat'),
            by = c('quantile', 'PFT')) %>% 
  arrange(PFT, quantile)


# functions ---------------------------------------------------------------


# functions for quantile matching using  fixed set of a few quantiles
qm_quant <- map(pfts, function(pft) {
  qm_quant_factory(from = quants$sw[[pft]],
                   to = quants$rap[[pft]])
  
})


# figures ----------------------------------------------------------------

# visual check that code is working (can compare to q-q plot made
# in 02_quantile_matching.R)

create_seq <- function(pft) {
  m <- max(sw1$biomass[sw1$PFT == pft])
  seq(from = 0, to = m + 10, length.out = 200)
}

df_seq <- tibble(Aherb = create_seq('Aherb'),
                 Pherb = create_seq('Pherb'))

seq1 <- df_seq %>% 
  mutate(obs = 1:nrow(.)) %>% 
  pivot_longer(-obs,
               names_to = "PFT",
               values_to = "biomass")
seq1$biomass_qm <- NA

for(pft in pfts) {
  seq1$biomass_qm[seq1$PFT == pft] <- qm_quant[[pft]](seq1$biomass[seq1$PFT == pft] )
}

ggplot(seq1, aes(biomass, biomass_qm)) +
  geom_line() +
  facet_wrap(~PFT, ncol = 1, scales = 'free') +
  geom_abline(color = 'gray') +
  geom_point(data = out, aes(biomass_stepwat, biomass_rap), color = 'blue')


# save objects ------------------------------------------------------------

run <- unique(sw1$run) 

stopifnot(length(run) == 1)

name <- paste0('quantiles_for_qm_', qual_cutoff, 'match', run, '.csv')

write_csv(out, file.path('data_processed/qm', name))
