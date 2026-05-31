#!/usr/bin/env Rscript
## Simulation_PriorSensitivity.R — Only Prior(II) needed (Prior(I) is the default used elsewhere)
## Prior(II) from paper:
##   A_k0=2, Lambda_k0=2, Sigma_Ak0=10^4*I_3, Sigma_Lk0=10^4*I_2,
##   a_eps0=3, b_eps0=2, beta_0=2, sigma_beta=10^4,
##   a_lambda1=4, a_lambda2=0.001
argv_file <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argv_file) == 0) argv_file <- "--file=Simulation_PriorSensitivity.R"
this_file <- sub("^--file=", "", argv_file[[1]])
source(file.path(dirname(this_file), "sim_common.R"))

defs <- list(
  list(
    tag = "prior_case_ii",
    label = "Prior(II), N=200, baseline constant, CR=30%",
    setting_overrides = list(
      n_subject = 200, baseline_mode = "constant", censoring_rate = 0.30,
      prior_variant = "II", missing_rate = 0, block_missing_prob = 0
    ),
    cli_overrides = list(
      ## Intercept A prior: A_k0 = 2, Sigma_Ak0 = 10^4 * I
      A_prior_mean = 2,
      A_prior_var = 10000,
      ## Loading prior: Lambda_k0 = 2, Sigma_Lk0 = 10^4 * I
      loading_prior_var = 10000,
      loading_prior_mean = 2,
      ## Psi prior: a_eps0 = 3, b_eps0 = 2
      psi_a0 = 3,
      psi_b0 = 2,
      ## Beta prior: beta_0 = 2, sigma_beta = 10^4 => sd = sqrt(10^4) = 100
      beta_prior_mean = 2,
      theta_prior_sd = 100,
      ## Lambda prior: a_lambda1 = 4, a_lambda2 = 0.001
      lambda_a1 = 4,
      lambda_a2 = 0.001
    )
  )
)

run_simulation_script("Simulation_PriorSensitivity", defs)
