#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_Main.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(tag = "main_const_N200_CR30", label = "constant baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_const_N500_CR30", label = "constant baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "constant", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_linear_N200_CR30", label = "linear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "linear", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_linear_N500_CR30", label = "linear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "linear", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_nonlinear_N200_CR30", label = "nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_nonlinear_N500_CR30", label = "nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_const_N200_CR50", label = "constant baseline, N=200, CR=50%", setting_overrides = list(n_subject = 200, baseline_mode = "constant", censoring_rate = 0.50, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "main_const_N500_CR50", label = "constant baseline, N=500, CR=50%", setting_overrides = list(n_subject = 500, baseline_mode = "constant", censoring_rate = 0.50, missing_rate = 0, block_missing_prob = 0))
)

run_simulation_script("Simulation_Main", defs)
