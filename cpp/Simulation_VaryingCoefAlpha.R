#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_VaryingCoefAlpha.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(tag = "varying_linear_N200", label = "alpha linear time-varying, N=200, baseline constant, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30, alpha_mode = "linear_tv", missing_rate = 0, block_missing_prob = 0)),
  list(tag = "varying_linear_N500", label = "alpha linear time-varying, N=500, baseline constant, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "constant", censoring_rate = 0.30, alpha_mode = "linear_tv", missing_rate = 0, block_missing_prob = 0)),
  list(tag = "varying_nonlinear_N200", label = "alpha nonlinear time-varying, N=200, baseline constant, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30, alpha_mode = "nonlinear_tv", missing_rate = 0, block_missing_prob = 0)),
  list(tag = "varying_nonlinear_N500", label = "alpha nonlinear time-varying, N=500, baseline constant, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "constant", censoring_rate = 0.30, alpha_mode = "nonlinear_tv", missing_rate = 0, block_missing_prob = 0))
)

run_simulation_script("Simulation_VaryingCoefAlpha", defs)
