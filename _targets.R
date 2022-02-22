library(targets)
# This is an example _targets.R file. Every
# {targets} pipeline needs one.
# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(summary) to view the results.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
source("functions_rpack.R")
source("CodeCollection/utility_functions.R")
source("CodeCollection/simulate_gamma_mixture.R")
source("CodeCollection/simulate_unif_grid.R")

# Set target-specific options such as packages.
tar_option_set(
  debug = "clust_dat4",
  packages = c(
    "tidyverse",
    "LaplacesDemon",
    "Matrix",
    "Gmedian",
    "lcmix",
    "flexmix",
    "purrr",
    "broman"
  )
)


# End this file with a list of target objects.
list(
  
  #### DATA SET ####
  tar_target(
    list_dat100,
    simulate_gamma_mixture_n(
      N = 10,
      n_total = 500,
      k = 10,
      n_out = 20,
      out_scale = 5,
      scale_between_range = c(0, 1),
      outgroup_alpha = 0.4,
      place_on_grid = TRUE,
      overlap_scale = 0.5
    )
  ),
  
  #### CLUSTERING PARAMETERS ####
  tar_target(mean_dat100, round(sum(list_dat100[[1]]$Y |> pull(w) / 10))),
  tar_target(range, c((mean_dat100 - 10*100), (mean_dat100 + 10*100))),
  tar_target(lambda_par, seq(0.1, 1, 0.1)),
  
  #### CLUSTERINGS ####
  tar_target(
    clust1,
    clust_with_params_list(
      dat_list = list_dat100,
      lambda_par = lambda_par,
      lambda = NULL,
      k = 10,
      N = 20,
      range = range,
      d = hav.dist2_par
    )
  ),
  tar_target(
    clust2,
    clust_with_params_list(
      dat_list = list_dat100,
      lambda_par = lambda_par,
      lambda = 0.01,
      k = 10,
      N = 20,
      range = range,
      d = hav.dist2_par
      
    )
  ),
  tar_target(
    clust3,
    clust_with_params_list(
      dat_list = list_dat100,
      lambda_par = lambda_par,
      lambda = 0.005,
      k = 10,
      N = 20,
      range = range,
      d = hav.dist2_par
    )
  ),
  tar_target(
    clust4,
    clust_with_params_list(
      dat_list = list_dat100,
      lambda_par = lambda_par,
      lambda = 0.001,
      k = 10,
      N = 20,
      range = range,
      d = hav.dist2_par
    )
  ),
  tar_target(
    clust5,
    clust_with_params_list(
      dat_list = list_dat100,
      lambda_par = lambda_par,
      lambda = 0.0001,
      k = 10,
      N = 20,
      range = range,
      d = hav.dist2_par
    )
  )
  
  
)
