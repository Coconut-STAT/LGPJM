#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_MainGPKernels.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(tag = "gp_main_const_N200_CR30", label = "GP omega truth, constant baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "constant", omega_mode = "gp_kernel", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_const_N500_CR30", label = "GP omega truth, constant baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "constant", omega_mode = "gp_kernel", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_linear_N200_CR30", label = "GP omega truth, linear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "linear", omega_mode = "gp_kernel", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_linear_N500_CR30", label = "GP omega truth, linear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "linear", omega_mode = "gp_kernel", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_nonlinear_N200_CR30", label = "GP omega truth, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", omega_mode = "gp_kernel", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_nonlinear_N500_CR30", label = "GP omega truth, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", omega_mode = "gp_kernel", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_const_N200_CR50", label = "GP omega truth, constant baseline, N=200, CR=50%", setting_overrides = list(n_subject = 200, baseline_mode = "constant", omega_mode = "gp_kernel", censoring_rate = 0.50, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "gp_main_const_N500_CR50", label = "GP omega truth, constant baseline, N=500, CR=50%", setting_overrides = list(n_subject = 500, baseline_mode = "constant", omega_mode = "gp_kernel", censoring_rate = 0.50, missing_rate = 0, block_missing_prob = 0))
)

run_simulation_script("Simulation_MainGPKernels", defs)
