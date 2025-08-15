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

# to run for all combinations of vr and years
scripts_vr_yr <- c(
  # 'scripts/07_bio_figs_weighted.R',
  'scripts/07_figs_sei_scd-adj.R'
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

# scripts that have years compontent (i.e. are run for only one time-period,
# at a time
# callr::rscript(scripts_vr_yr, cmdargs = uu(arg_list_end)) 
#callr::rscript(scripts_vr_yr, cmdargs = uu(arg_list_mid))



# have run
# callr::rscript('scripts/05_interpolated_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_end)) 
# callr::rscript('scripts/05_interpolated_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_mid)) 
# callr::rscript('scripts/06_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_mid))
# callr::rscript('scripts/07_figs_sei_scd-adj.R', cmdargs = uu(arg_list_mid)) # done
# callr::rscript('scripts/07_figs_sei_scd-adj.R', cmdargs = uu(arg_list_end)) 

for(script in scripts_vr_yr) {
  callr::rscript(script, cmdargs = uu(arg_list_end))
  callr::rscript(script, cmdargs = uu(arg_list_mid))
  callr::rscript(script, cmdargs = uu(arg_list_r1.0_end))
  callr::rscript(script, cmdargs = uu(arg_list_r1.0_mid))
}

# r1.0 --------------------------------------------------------------------

arg_list_mid$vr <- "--vr=r1.0"
# callr::rscript('scripts/05_interpolated_summarize_sei_scd-adj.R', 
#                cmdargs = uu(arg_list_r1.0_mid)) 
# callr::rscript('scripts/06_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_mid))
# callr::rscript('scripts/06_summarize_sei_scd-adj.R', cmdargs = uu(arg_list_mid))
#callr::rscript('scripts/07_figs_sei_scd-adj.R', cmdargs = uu(arg_list_mid)) 