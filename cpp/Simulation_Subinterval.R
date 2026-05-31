#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_Subinterval.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(tag = "subinterval_G5_N200", label = "G=5, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, baseline_partition_G = 5, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "subinterval_G5_N500", label = "G=5, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, baseline_partition_G = 5, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "subinterval_G10_N200", label = "G=10, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, baseline_partition_G = 10, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "subinterval_G10_N500", label = "G=10, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, baseline_partition_G = 10, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "subinterval_G15_N200", label = "G=15, nonlinear baseline, N=200, CR=30%", setting_overrides = list(n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30, baseline_partition_G = 15, missing_rate = 0, block_missing_prob = 0)),
  list(tag = "subinterval_G15_N500", label = "G=15, nonlinear baseline, N=500, CR=30%", setting_overrides = list(n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30, baseline_partition_G = 15, missing_rate = 0, block_missing_prob = 0))
)

run_simulation_script("Simulation_Subinterval", defs)
