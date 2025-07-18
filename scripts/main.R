# for running all scripts 
# including some scripts repeatedly with different parameters

# (e.g. re-rerun for different 'runs', time-periods, or ecoregions)
# see src/params for defaults


# list scripts ------------------------------------------------------------

# scripts that need to be run for different time-periods
scripts_vr_yr <- c(
  'scripts/05_interpolated_summarize_sei_scd-adj.R'
  # 'scripts/06_summarize_sei_scd-adj.R',
)

# scripts the need to be run seperately for each version of the
# ecoregions, # but not years
scripts_vr <- c(
  # 'scripts/04_interpolation_weights.R'
  # 'scripts/04_summarize_fire_drivers.R'
  # 'scripts/06_fire_area.R',
)

# parameters --------------------------------------------------------------

end <- "--years=2070-2100"
mid <- "--years=2030-2060"
arg_list_end <- list(
  vr = "--vr=r1.1", # ecoregion version
  years = end
)
arg_list_mid <- arg_list_end
arg_list_mid$years <- mid

# functions ---------------------------------------------------------------

uu <- function(x) unname(unlist(x))

# run scripts -------------------------------------------------------------

# scripts that have years compontent (i.e. are run for only one time-period,
# at a time
# callr::rscript(scripts_vr_yr, cmdargs = uu(arg_list_end)) 
#callr::rscript(scripts_vr_yr, cmdargs = uu(arg_list_mid))

# callr::rscript('scripts/05_interpolated_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_end))

# still need to run
# callr::rscript('scripts/05_interpolated_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_mid))
# callr::rscript('scripts/06_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_mid))
