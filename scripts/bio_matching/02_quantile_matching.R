# Purpose: fit ecdfs to RAP and stepwat biomass so can do quantile matching

# Martin Holdrege

# Started: Nov 8, 2024

# params ------------------------------------------------------------------

# quantile matching functions fit using two dataset and two methods
# but only one is outputted, specified here
dataset_abbrev <- 'interp'
method_abbrev <- 'ecdf'

# dependencies ------------------------------------------------------------

library(tidyverse)
source('src/fig_functions.R')
source('src/qm_functions.R')
theme_set(theme_custom1())

# read in data ------------------------------------------------------------

# list created in 01_combine_data.R
m <- readRDS('data_processed/temp_rds/rap_sw_matching.rds') 
run <- "fire1_eind1_c4grass1_co20_2311"

# vectors -----------------------------------------------------------------

pfts <- c('Pherb' = 'Pherb', 'Aherb' = 'Aherb')
stopifnot(run %in% m$runs)

# fitting cdfs ------------------------------------------------------------

# rap datasets
rap_df <- m$rap_comb

sw_df <- m$sw_comb %>% 
  filter(run == !!run)

# fit emperical cdfs to interpolated data
cdfs_interp <- list()
dataset <- 'interpolated'
cdfs_interp$rap <- map(pfts, function(pft) {
  x <- rap_df[rap_df$PFT == pft  & rap_df$dataset == dataset, ]$biomass
  ecdf(x)
})

cdfs_interp$sw <- map(pfts, function(pft) {
  x <- sw_df[sw_df$PFT == pft & sw_df$dataset == dataset, ]$biomass
  ecdf(x)
})

# ecdfs fit to site level (not interpolated data)
cdfs_site <- list()
cdfs_site$rap <- map(pfts, function(pft) {
  x <- rap_df[rap_df$PFT == pft & rap_df$dataset != dataset, ]$biomass
  ecdf(x)
})

cdfs_site$sw <- map(pfts, function(pft) {
  x <- sw_df[sw_df$PFT == pft & sw_df$dataset != dataset, ]$biomass
  ecdf(x)
})


# functions ---------------------------------------------------------------

qm_factory <- function(pft, cdfs) {
  # using this approach so the ecdfs are encapsulated in the function environment
  cdf_from <- cdfs$sw[[pft]] 
  cdf_to <- cdfs$rap[[pft]]
  function(x) {
    probs <- cdf_from(x)
    values <- quantile(cdf_to, probs = probs)
    as.numeric(values)
  }
}

# quantile matching from fixed set of a few quantiles
qm_quant_factory <- function(from, to) {
  function(x) {
    qm_from_quantiles(x, from = from, to = to)
  }
}

probs = c(seq(0, 0.9, by = 0.1),
          0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1)# probabilities to calculate the quantiles of

# quantiles
quants_interp <- map(cdfs_interp, function(l) {
  map(l, \(x) as.numeric(quantile(x, probs = probs)))
})

quants_site <-map(cdfs_site, function(l) {
  map(l, \(x) as.numeric(quantile(x, probs = probs)))
})


# create functions where input is a vector of biomass
# and output is a vector of adjusted biomass, so that it
# is on the same scale as RAP
qm_interp <- map(pfts, \(x) qm_factory(x, cdfs_interp)) # list of two functions
qm_site <- map(pfts, \(x) qm_factory(x, cdfs_site))

# functions for quantile matching using  fixed set of a few quantiles
qm_quant_interp <- map(pfts, function(pft) {
  qm_quant_factory(from = quants_interp$sw[[pft]],
                   to = quants_interp$rap[[pft]])
  
})

qm_quant_site <- map(pfts, function(pft) {
  qm_quant_factory(from = quants_site$sw[[pft]],
                   to = quants_site$rap[[pft]])
  
})

# list of list of functions
qm_comb <- list(
  'ecdf_interp' = qm_interp,
  'ecdf_site' = qm_site,
  'quant_interp' = qm_quant_interp,
  'quant_site' = qm_quant_site
)

# figures ----------------------------------------------------------------

# sequences of original and quantile matched biomass values, for plotting
create_seq <- function(pft) {
  m <- max(sw_df$biomass[sw_df$PFT == pft])
  seq(from = 0, to = m + 10, length.out = 200)
}

df_seq <- tibble(Aherb = create_seq('Aherb'),
                 Pherb = create_seq('Pherb'))

dat <- map2_dfc(names(qm_comb), qm_comb, function(name, l) {
  dat <- map_dfc(pfts, function(pft) {
    f <- l[[pft]]
    f(df_seq[[pft]])
  })
  names(dat) <- paste0(pfts, "_", name)
  dat
})

obs <- 1:nrow(df_seq)
seq_long <- df_seq %>% 
  mutate(obs = obs) %>% 
  pivot_longer(-obs,
               names_to = "PFT",
               values_to = "biomass")

lookup_method <- c(ecdf = 'ecdf function',
  quant = 'fixed quantiles')

lookup_dataset <- c(interp = 'interpolated',
                    site = 'site level')

df_seq2 <- dat %>% 
  mutate(obs = obs) %>% 
  pivot_longer(
    cols = matches("(.*)_(.*)_(.*)"),
    names_to = c("PFT", "method", "dataset"),
    names_pattern = "(.*)_(.*)_(.*)",
    values_to = 'biomass_qm'
  ) %>% 
  left_join(seq_long, by = c("PFT", "obs")) %>% 
  select(-obs) %>% 
  # using lookup vectors to change values
  mutate(dataset = lookup_dataset[dataset],
         method = lookup_method[method])


pdf("figures/bio_matching/q-q_plots_v1.pdf")
ggplot(df_seq2, aes(x = biomass, y = biomass_qm)) +
  geom_abline(slope = 1, color = 'gray') +
  geom_line(aes(linetype = method, color = dataset)) +
  facet_wrap(~PFT, scales = 'free', ncol = 1) +
  labs(x = 'Biomass (stepwat scale)',
       y = 'Quantile matched to RAP scale',
       caption = paste0(
         length(probs), ' fixed quantiles used.',
         "\n all grazing levels from ", run, " used"
       )) 
dev.off()

# save objects ------------------------------------------------------------

name <- paste0(method_abbrev, "_", dataset_abbrev)

dataset_descript <- if(dataset_abbrev == "interp"){
  paste0(lookup_dataset[dataset_abbrev],' (<', m$qual_cutoff, 
         'matching quality)')
} else {
  lookup_dataset[dataset_abbrev]
}

out <- list('Aherb' = qm_comb[[name]]$Aherb,
            'Pherb' = qm_comb[[name]]$Pherb,
            # update this caption depending on what qm functions are put here
            caption = paste(lookup_method[method_abbrev], 
                            'used for quantile mapping,',
                            '\n with',
                            dataset_descript,
                            'data from', run),
            'dataset' = lookup_dataset[dataset_abbrev],
            qual_cutoff = m$qual_cutoff,
            run = run)
saveRDS(out, 'data_processed/temp_rds/qm_funs.rds')
