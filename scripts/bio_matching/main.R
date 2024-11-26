# run other scripts

# List of scripts to source
scripts <- c('scripts/bio_matching/01_combine_data.R', 
             'scripts/bio_matching/02_quantile_matching.R', 
             'scripts/bio_matching/03_distribution_figures.R', 
             'scripts/bio_matching/03_fire_derivative.R', 
             'scripts/bio_matching/03_predict_fire-prob.R')

# Loop through each script and source it in its own environment
for (script in scripts) {
  print(script)
  
  # Create a new isolated environment
  # script_env <- new.env()
  
  # Source the script in the new environment
  #sys.source(file = script, envir = script_env, toplevel.env = script_env)
  source(script)
  # (Optional) Print the environment content for debugging
  # print(ls(envir = script_env))
  rm(list = ls())
  # Remove the temporary environment to free up memory
  #rm(script_env)
  
  # Call garbage collection to release memory
  #gc()
}

