#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_Missing.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(tag = "missing10_N200", label = "missing rate 10%, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_mechanism = "MAR", missing_rate = 0.10, block_missing_prob = 0.03)),
  list(tag = "missing10_N500", label = "missing rate 10%, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_mechanism = "MAR", missing_rate = 0.10, block_missing_prob = 0.03)),
  list(tag = "missing30_N200", label = "missing rate 30%, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_mechanism = "MAR", missing_rate = 0.30, block_missing_prob = 0.03)),
  list(tag = "missing30_N500", label = "missing rate 30%, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_mechanism = "MAR", missing_rate = 0.30, block_missing_prob = 0.03)),
  list(tag = "missing50_N200", label = "missing rate 50%, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_mechanism = "MAR", missing_rate = 0.50, block_missing_prob = 0.03)),
  list(tag = "missing50_N500", label = "missing rate 50%, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, missing_mechanism = "MAR", missing_rate = 0.50, block_missing_prob = 0.03))
)

run_simulation_script("Simulation_Missing", defs)
