#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_ModelComparison.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(
    tag = "model_case_a_gp",
    label = "Case (a): linear omega truth, GPJM model",
    methods = "internal",
    cli_overrides = list(compute_criteria = "true"),
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      omega_mode = "linear", missing_rate = 0, block_missing_prob = 0
    )
  ),
  list(
    tag = "model_case_a_linear",
    label = "Case (a): linear omega truth, linear model",
    methods = "linear_model",
    cli_overrides = list(compute_criteria = "true"),
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      omega_mode = "linear", missing_rate = 0, block_missing_prob = 0
    )
  ),
  list(
    tag = "model_case_b_gp",
    label = "Case (b): nonlinear omega truth, GPJM model",
    methods = "internal",
    cli_overrides = list(compute_criteria = "true"),
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      omega_mode = "sinusoidal", missing_rate = 0, block_missing_prob = 0
    )
  ),
  list(
    tag = "model_case_b_linear",
    label = "Case (b): nonlinear omega truth, linear model",
    methods = "linear_model",
    cli_overrides = list(compute_criteria = "true"),
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      omega_mode = "sinusoidal", missing_rate = 0, block_missing_prob = 0
    )
  )
)

run_simulation_script("Simulation_ModelComparison", defs)
