# main parameters used in other scripts
# (so they don't need to be defined in each script)


# options can be passed to this script, or scripts sourcing it
# to change these defaults
# callr::rscript('scripts/nameofscript.R', cmdargs = c('--vr=r1.0', 
#    '--years=2030-2060')
option_list <- list(
  optparse::make_option("--run", type = "character", 
              default = "fire1_eind1_c4grass1_co20_2503"),
  # interpolation version
  optparse::make_option("--v_interp", type = "character", default = "v4"),
  # version of region polygons to use (summaries made by regions)
  # r1.0 is 4 regions
  # r1.1 is 8 regions
  optparse::make_option("--vr", type = "character", default = "r1.0"),
  optparse::make_option("--years", type = "character", default = "2030-2060"),
  # included --test_run with no value forst test_run to be true
  # otherwise leave it out
  optparse::make_option("--test_run", action = "store_true", default = FALSE)
)

opt_parser <-optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
run <- opt$run

v_interp <- opt$v_interp

vr <- opt$vr

years <- opt$years

stopifnot(years %in% c('2030-2060', '2070-2100'))

# for file name
opt$vr_name <- if(vr == 'r1.0') "" else paste0("_", vr)
vr_name <- opt$vr_name

if(years == '2070-2100') {
  opt$yr_lab <- '' # previously end was done by default so not changing that in the filename
} else if (years == '2030-2060') {
  opt$yr_lab <- '_mid'
}

print(opt)