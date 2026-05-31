#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_ErrorSensitivity.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(
    tag = "error_case_i_t3_scaled",
    label = "Case (i): error = scaled t(3), N=200, baseline constant, CR=30%",
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      error_mode = "t3_scaled", missing_rate = 0, block_missing_prob = 0
    )
  ),
  list(
    tag = "error_case_ii_shifted_gamma",
    label = "Case (ii): error = shifted gamma, N=200, baseline constant, CR=30%",
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      error_mode = "shifted_gamma", missing_rate = 0, block_missing_prob = 0
    )
  )
)

run_simulation_script("Simulation_ErrorSensitivity", defs)
