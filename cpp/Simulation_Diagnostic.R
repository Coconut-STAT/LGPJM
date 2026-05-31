#!/usr/bin/env Rscript
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_Diagnostic.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(
    tag = "diag_main_nonlinear_N500_CR30",
    label = "diagnostic setting: nonlinear baseline, N=500, CR=30%, 3 chains, 1 replicate",
    chains = 3,
    datasets = 1,
    setting_overrides = list(
      n_subject = 500, baseline_mode = "nonlinear", censoring_rate = 0.30,
      omega_mode = "sinusoidal", missing_rate = 0, block_missing_prob = 0
    ),
    cli_overrides = list(save_chains = "true")
  )
)

run_simulation_script("Simulation_Diagnostic", defs)
