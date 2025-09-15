# for running all scripts 
# including some scripts repeatedly with different parameters

# (e.g. re-rerun for different 'runs', time-periods, or ecoregions)
# see src/params for defaults


# list scripts ------------------------------------------------------------


# scripts the need to be run seperately for each version of the
# ecoregions, # but not years
scripts_vr <- c(
  # 'scripts/04_interpolation_weights.R'
  # 'scripts/04_summarize_fire_drivers.R'
  #'scripts/06_fire_area.R',
  # 'scripts/07_fire_area_summaries.R'
)

# scripts that need to be run seperately for each set of years
# but not ecoregion versions
scripts_yr <- c(
  "scripts/06_sei_attribution.R"
)

# to run for all combinations of vr and years
scripts_vr_yr <- c(
  #'scripts/05_interpolated_summarize_sei_scd-adj.R'
  # 'scripts/06_summarize_sei_scd-adj.R',
  # 'scripts/07_bio_figs_weighted.R',
  #'scripts/07_figs_sei_scd-adj.R'
  #'scripts/08_fire_area_figs.R'
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
arg_list_r1.0_mid <- arg_list_mid
arg_list_r1.0_mid$vr <- "--vr=r1.0"
arg_list_r1.0_end <- arg_list_end
arg_list_r1.0_end$vr <- "--vr=r1.0"
# functions ---------------------------------------------------------------

uu <- function(x) unname(unlist(x))

# run scripts -------------------------------------------------------------


if(length(scripts_vr) > 0) {
  for(script in scripts_vr) {
    callr::rscript(script, cmdargs = uu(arg_list_r1.0_end))
    callr::rscript(script, cmdargs = uu(arg_list_end))
  }
}

if(length(scripts_yr) > 0) {
  for(script in scripts_yr) {
    callr::rscript(script, cmdargs = uu(arg_list_end))
    callr::rscript(script, cmdargs = uu(arg_list_mid))
  }
}

if(length(scripts_vr_yr) > 0) {
  for(script in scripts_vr_yr) {
    callr::rscript(script, cmdargs = uu(arg_list_end))
    callr::rscript(script, cmdargs = uu(arg_list_mid))
    callr::rscript(script, cmdargs = uu(arg_list_r1.0_end))
    callr::rscript(script, cmdargs = uu(arg_list_r1.0_mid))
  }
}
