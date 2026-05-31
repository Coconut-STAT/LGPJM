get_default_loading <- function() {
  mat <- matrix(0, nrow = 9, ncol = 2)
  mat[1:5, 1] <- 0.8
  mat[6, 1] <- 1
  mat[7, 2] <- 1
  mat[8:9, 2] <- 0.8
  mat
}

get_alpha_function <- function(alpha_mode) {
  switch(
    alpha_mode,
    constant = function(t) {
      rbind(rep(0.5, length(t)), rep(0.5, length(t)))
    },
    linear_tv = function(t) {
      rbind(-0.3 * t + 0.8, 0.3 * t + 0.2)
    },
    nonlinear_tv = function(t) {
      rbind(0.25 * sin(2 * t) + 0.5, 0.25 * cos(2 * t) + 0.5)
    },
    stop(sprintf("Unsupported alpha mode: %s", alpha_mode))
  )
}

get_baseline_hazard <- function(baseline_mode) {
  switch(
    baseline_mode,
    constant = function(t) rep(1, length(t)),
    linear = function(t) t + 1,
    nonlinear = function(t) 0.5 * t^2 + 1,
    stop(sprintf("Unsupported baseline mode: %s", baseline_mode))
  )
}

apply_setting_functions <- function(st) {
  st$alpha_fun <- get_alpha_function(st$alpha_mode)
  st$baseline_fun <- get_baseline_hazard(st$baseline_mode)
  st$time_points <- seq(st$t_start, st$t_end, length.out = st$n_time)
  st
}

get_common_setting <- function() {
  list(
    n_subject = 200,
    n_time = 10,
    t_start = 0,
    t_end = 2,
    p = 9,
    q = 2,
    h = 3,
    ordinal_categories = c(2, 3, 4),
    loading_true = get_default_loading(),
    loading_group = c(rep(1, 6), rep(2, 3)),
    fixed_loading_index = c(6, 7),
    A_true = matrix(rep(c(-2, -2, 1), 9), nrow = 9, byrow = TRUE),
    psi_true = c(rep(1, 3), rep(0.3, 6)),
    beta_true = c(0.5, -0.5, 0.5),
    alpha_mode = "constant",
    baseline_mode = "constant",
    omega_mode = "sinusoidal",
    omega_kernel = list(sigma2 = 1, length_scale = 0.5, nugget = 1e-6),
    censoring_rate = 0.30,
    missing_mechanism = "MAR",
    missing_rate = 0.30,
    block_missing_prob = 0.03,
    error_mode = "normal",
    prior_variant = "I",
    baseline_partition_method = "quantile",
    baseline_partition_G = 5,
    ## ----- New MCMC parameters (matching reference) -----
    N_int = 20,                    # integration grid size for survival integral
    sigmaf_omega = 1,              # GP kernel sigmaf for omega interpolation
    l_omega = 2,                   # GP kernel length scale for omega interpolation
    l_B = 0.25,                    # GP kernel length scale for B matrix (Qnew)
    sigma_alpha_ph = 1,            # GP kernel sigmaf for alpha_ph prior
    l_alpha_ph = 1,                # GP kernel length scale for alpha_ph prior
    gp_pred_nugget = 0.1,          # nugget used in GP interpolation to integration grid
    c_omega = 0.5,                 # omega MH proposal scale
    c_omega_linear = 1.0,          # linear omega independence-proposal scale
    c_alpha = 1.0,                 # alpha_ph MH proposal scale for stable GP block proposal
    c_beta = 2,                    # beta MH proposal scale (reference default)
    sigma_tau = 0.1,               # threshold proposal base variance; fit config divides by n_time
    l_B_lower = 0.02,              # lower bound for latent GP length-scale MH
    l_B_upper = 5,                 # upper bound for latent GP length-scale MH
    l_B_proposal_width = 0.20      # multiplicative U(1-a_l, 1+a_l) proposal width
  )
}
