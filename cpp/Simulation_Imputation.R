#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_Imputation.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(
    tag = "imputation_internal_N200",
    label = "internal handling, nonlinear baseline, N=200, CR=30%",
    methods = "internal",
    setting_overrides = list(
      n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30,
      missing_mechanism = "MAR", missing_rate = 0.30, block_missing_prob = 0.03
    )
  ),
  list(
    tag = "imputation_locf_N200",
    label = "LOCF imputation, nonlinear baseline, N=200, CR=30%",
    methods = "imputation",
    setting_overrides = list(
      n_subject = 200, baseline_mode = "nonlinear", censoring_rate = 0.30,
      missing_mechanism = "MAR", missing_rate = 0.30, block_missing_prob = 0.03
    )
  ),
  list(
    tag = "imputation_internal_N500",
    label = "internal handling, nonlinear baseline, N=500, CR=30%",
    methods = "internal",
    setting_overrides = list(
      n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30,
      missing_mechanism = "MAR", missing_rate = 0.30, block_missing_prob = 0.03
    )
  ),
  list(
    tag = "imputation_locf_N500",
    label = "LOCF imputation, nonlinear baseline, N=500, CR=30%",
    methods = "imputation",
    setting_overrides = list(
      n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30,
      missing_mechanism = "MAR", missing_rate = 0.30, block_missing_prob = 0.03
    )
  )
)

run_simulation_script("Simulation_Imputation", defs)
