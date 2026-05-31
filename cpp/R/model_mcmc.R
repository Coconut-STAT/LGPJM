make_baseline_cuts <- function(obs_time, delta, G, method, t_start, t_end) {
  method <- tolower(method)
  t_end_use <- max(t_end, obs_time, na.rm = TRUE)
  if (method == "equal") return(seq(t_start, t_end_use, length.out = G + 1))
  event_t <- obs_time[delta == 1]
  if (length(event_t) < G + 1) return(seq(t_start, t_end_use, length.out = G + 1))
  qv <- as.numeric(stats::quantile(event_t, probs = seq(0, 1, length.out = G + 1), names = FALSE))
  qv[1] <- t_start
  qv[length(qv)] <- t_end_use
  for (g in 2:length(qv)) if (qv[g] <= qv[g - 1]) qv[g] <- qv[g - 1] + 1e-4
  qv
}

require_lgpjm_cpp <- function(name) {
  if (!exists(name, mode = "function")) {
    stop(sprintf("Required C++ sampler function `%s` is not loaded. Re-run Rcpp::sourceCpp on cpp/cpp/mcmc_core.cpp.", name))
  }
}

init_thresholds <- function(data, setting) {
  n_ord <- length(setting$ordinal_categories)
  out <- vector("list", n_ord)
  for (j in seq_len(n_ord)) {
    K <- setting$ordinal_categories[j]
    if (!is.null(data$thresholds) && length(data$thresholds) >= j && length(data$thresholds[[j]]) == K - 1) {
      tau <- as.numeric(data$thresholds[[j]])
    } else {
      obs <- as.numeric(data$z_ord_obs[j, , ])
      obs <- obs[!is.na(obs)]
      if (length(obs) == 0) {
        probs <- seq(0, 1, length.out = K + 1)[2:K]
      } else {
        probs <- vapply(seq_len(K - 1), function(k) mean(obs <= k), numeric(1))
      }
      tau <- stats::qnorm(pmin(pmax(probs, 1e-4), 1 - 1e-4))
    }
    if (length(tau) > 1) {
      for (r in 2:length(tau)) if (tau[r] <= tau[r - 1]) tau[r] <- tau[r - 1] + 0.25
    }
    out[[j]] <- tau
  }
  out
}

threshold_param_names <- function(thresholds) {
  names <- character()
  for (j in seq_along(thresholds)) {
    tau <- thresholds[[j]]
    if (length(tau) > 1) names <- c(names, paste0("tau", j, "_", 2:length(tau)))
  }
  names
}

flatten_free_thresholds <- function(thresholds) {
  vals <- numeric()
  for (tau in thresholds) if (length(tau) > 1) vals <- c(vals, tau[2:length(tau)])
  vals
}

thresholds_to_matrix <- function(thresholds) {
  lens <- as.integer(lengths(thresholds))
  max_len <- max(lens, 1L)
  mat <- matrix(0, nrow = length(thresholds), ncol = max_len)
  for (j in seq_along(thresholds)) {
    if (lens[j] > 0L) mat[j, seq_len(lens[j])] <- as.numeric(thresholds[[j]])
  }
  list(mat = mat, len = lens)
}

thresholds_from_matrix <- function(mat, lens) {
  lapply(seq_along(lens), function(j) {
    if (lens[j] <= 0L) return(numeric())
    as.numeric(mat[j, seq_len(lens[j])])
  })
}

build_Qnew_cube_from_l <- function(setting, l_vec) {
  Tn <- length(setting$time_points)
  q <- length(l_vec)
  out <- array(0, dim = c(Tn + 1, Tn + 1, q))
  for (k in seq_len(q)) {
    Q <- se_kernel_mat(setting$time_points, setting$time_points,
                       setting$sigmaf_omega, l_vec[k])
    out[, , k] <- build_Qnew(Q, Tn) + diag(setting$omega_kernel$nugget %||% 1e-6, Tn + 1)
  }
  out
}

invert_Qnew_cube <- function(Qnew) {
  out <- Qnew
  for (k in seq_len(dim(Qnew)[3])) {
    out[, , k] <- tryCatch(solve(Qnew[, , k]), error = function(e) MASS::ginv(Qnew[, , k]))
  }
  out
}

compute_kesi_p_state <- function(omega, data, setting, fit_cfg, l_vec = NULL) {
  if (identical(fit_cfg$fit_model, "linear")) {
    if (exists("cpp_lgpjm_compute_kesi_p_linear", mode = "function")) {
      return(cpp_lgpjm_compute_kesi_p_linear(omega, data$obs_time, setting$time_points,
        as.integer(setting$N_int), setting$t_start))
    }
    return(compute_all_kesi_p_linear(omega, data$obs_time, setting$time_points,
      setting$N_int, setting$t_start))
  }
  if (exists("cpp_lgpjm_compute_kesi_p", mode = "function")) {
    if (is.null(l_vec)) l_vec <- rep(setting$l_omega, setting$q)
    return(cpp_lgpjm_compute_kesi_p(omega, data$obs_time, setting$time_points,
      as.integer(setting$N_int), setting$t_start, setting$sigmaf_omega,
      as.numeric(l_vec), fit_cfg$gp_pred_nugget))
  }
  if (is.null(l_vec)) l_vec <- rep(setting$l_omega, setting$q)
  compute_all_kesi_p(omega, data$obs_time, setting$time_points, setting$N_int,
                     setting$t_start, setting$sigmaf_omega, l_vec[1])
}

compute_alpha_ph_p_state <- function(alpha_ph, data, setting, fit_cfg) {
  if (exists("cpp_lgpjm_compute_alpha_ph_p", mode = "function")) {
    return(cpp_lgpjm_compute_alpha_ph_p(alpha_ph, data$obs_time, setting$time_points,
      as.integer(setting$N_int), setting$t_start, setting$sigma_alpha_ph,
      setting$l_alpha_ph, fit_cfg$gp_pred_nugget))
  }
  compute_all_alpha_ph_p(alpha_ph, data$obs_time, setting$time_points, setting$N_int,
                         setting$t_start, setting$sigma_alpha_ph, setting$l_alpha_ph)
}

fit_linear_omega_state <- function(omega_obs, obs_time, time_points) {
  q <- dim(omega_obs)[1]
  Tn <- min(dim(omega_obs)[2], length(time_points))
  n <- dim(omega_obs)[3]
  t_use <- time_points[seq_len(Tn)]
  X <- cbind(1, t_use)
  XtX <- crossprod(X) + 1e-8 * diag(2)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
  omega_full <- array(0, dim = c(q, Tn + 1, n))
  omega_full[, seq_len(Tn), ] <- omega_obs[, seq_len(Tn), , drop = FALSE]
  for (i in seq_len(n)) {
    for (k in seq_len(q)) {
      y <- as.numeric(omega_obs[k, , i])
      coef <- XtX_inv %*% crossprod(X, y)
      omega_full[k, Tn + 1, i] <- coef[1] + coef[2] * obs_time[i]
    }
  }
  omega_full
}

init_chain_state <- function(data, setting, fit_cfg, chain_seed) {
  set.seed(chain_seed)
  n <- dim(data$z_ord_obs)[3]
  Tobs <- length(setting$time_points)
  Ty <- dim(data$z_ord_obs)[2]
  Tall <- max(Tobs + 1L, Ty)
  n_ord <- dim(data$z_ord_obs)[1]
  n_cont <- dim(data$x_cont_obs)[1]
  p <- setting$p
  n_cov <- ncol(data$d_cov)
  thresholds <- init_thresholds(data, setting)

  y_star <- array(NA_real_, dim = c(p, Ty, n))
  for (j in seq_len(n_ord)) {
    bounds <- c(-Inf, thresholds[[j]], Inf)
    for (i in seq_len(n)) for (m in seq_len(Ty)) {
      z <- data$z_ord_obs[j, m, i]
      if (is.na(z)) {
        y_star[j, m, i] <- rnorm(1)
      } else {
        lo <- bounds[z]
        up <- bounds[z + 1]
        y_star[j, m, i] <- if (is.finite(lo) && is.finite(up)) {
          0.5 * (lo + up)
        } else if (!is.finite(lo)) {
          up - 1
        } else {
          lo + 1
        }
      }
    }
  }
  for (j in seq_len(n_cont)) {
    gj <- n_ord + j
    for (i in seq_len(n)) for (m in seq_len(Ty)) {
      x <- data$x_cont_obs[j, m, i]
      y_star[gj, m, i] <- if (is.na(x)) rnorm(1) else x
    }
  }

  loading <- numeric(p)
  loading[setting$fixed_loading_index] <- 1
  free_idx <- setdiff(seq_len(p), setting$fixed_loading_index)
  loading[free_idx] <- rnorm(length(free_idx), mean = 0.5, sd = 0.2)
  psi <- c(rep(1, n_ord), rep(0.5, n_cont))
  omega <- array(rnorm(setting$q * Tall * n, sd = 0.25), dim = c(setting$q, Tall, n))
  beta <- rnorm(setting$h, sd = 0.15)
  alpha_ph <- matrix(0, nrow = setting$q, ncol = Tobs)
  lambda <- rep(0.2, fit_cfg$G)

  A <- setting$A_true
  for (j in seq_len(p)) {
    y_vec <- y_star[j, 1, ]
    fit <- tryCatch(lm.fit(data$d_cov[, -1, drop = FALSE], y_vec - A[j, 1]), error = function(e) NULL)
    if (!is.null(fit) && !any(is.na(fit$coefficients))) A[j, 2:n_cov] <- fit$coefficients
  }
  A[, 1] <- setting$A_true[, 1]

  if (identical(fit_cfg$fit_model, "linear")) {
    if (exists("cpp_update_omega_linear", mode = "function")) {
      omega_obs <- array(0, dim = c(setting$q, Tobs, n))
      omega_obs <- tryCatch(
        cpp_update_omega_linear(omega_obs, y_star, loading, psi, A, setting$time_points,
          as.integer(setting$loading_group), setting$q, fit_cfg$omega_linear_prior_var),
        error = function(e) NULL
      )
      if (is.null(omega_obs)) omega_obs <- omega[, seq_len(Tobs), , drop = FALSE]
    } else {
      omega_obs <- omega[, seq_len(Tobs), , drop = FALSE]
    }
    omega <- fit_linear_omega_state(omega_obs, data$obs_time, setting$time_points)
  }

  l_omega <- rep(setting$l_B, setting$q)
  Qnew <- build_Qnew_cube_from_l(setting, l_omega)
  Qnew_inv <- invert_Qnew_cube(Qnew)
  theta_phi <- rep(1, setting$q)
  B_diag <- rep(1, setting$q)
  kesi_p <- compute_kesi_p_state(omega, data, setting, fit_cfg, rep(setting$l_omega, setting$q))
  alpha_ph_p <- compute_alpha_ph_p_state(alpha_ph, data, setting, fit_cfg)

  list(
    y_star = y_star, loading = loading, psi = psi, omega = omega,
    beta = beta, alpha_ph = alpha_ph, lambda = lambda, A = A,
    thresholds = thresholds, Qnew = Qnew, Qnew_inv = Qnew_inv,
    theta_phi = theta_phi, B_diag = B_diag, l_omega = l_omega, kesi_p = kesi_p,
    alpha_ph_p = alpha_ph_p, accept_omega = 0L, total_omega = 0L,
    accept_alpha = 0L, total_alpha = 0L, accept_beta = 0L, total_beta = 0L,
    accept_threshold = 0L, total_threshold = 0L, accept_l = 0L, total_l = 0L
  )
}

compute_A_d <- function(state, data) state$A %*% t(data$d_cov)

omega_regular_matrix <- function(state, k, Tn) {
  matrix(state$omega[k, seq_len(Tn), ], nrow = Tn, ncol = dim(state$omega)[3])
}

update_gp_hyperparams <- function(state, data, setting, fit_cfg) {
  require_lgpjm_cpp("cpp_lgpjm_gp_hyper_update")
  out <- cpp_lgpjm_gp_hyper_update(state$omega, state$kesi_p, state$Qnew_inv,
    state$B_diag, state$theta_phi, state$l_omega, state$alpha_ph_p, state$beta,
    state$lambda, data$obs_time, as.integer(data$delta), data$u_cov, state$cuts,
    setting$time_points, as.integer(setting$N_int), setting$t_start,
    setting$sigmaf_omega, setting$omega_kernel$nugget %||% 1e-6,
    fit_cfg$gp_sigma2_a0, fit_cfg$gp_sigma2_b0, fit_cfg$gp_l_prop_width,
    fit_cfg$gp_l_lower, fit_cfg$gp_l_upper, state$accept_l, state$total_l)
  state$theta_phi <- out$theta_phi
  state$B_diag <- out$B_diag
  state$l_omega <- out$l_vec
  state$Qnew_inv <- out$Qnew_inv
  state$kesi_p <- out$kesi_p
  state$accept_l <- out$accept_l
  state$total_l <- out$total_l
  state
}

update_A <- function(state, data, setting, fit_cfg) {
  require_lgpjm_cpp("cpp_lgpjm_update_A_scaled")
  state$A <- cpp_lgpjm_update_A_scaled(state$A, state$y_star, state$omega, state$loading,
    state$psi, data$d_cov, as.integer(setting$loading_group), fit_cfg$A_prior_mean,
    fit_cfg$Sigma_A0_inv, setting$A_true[, 1])
  state
}

update_loadings <- function(state, data, setting, prior_var = 100, prior_mean = 0, A_d = NULL) {
  require_lgpjm_cpp("cpp_lgpjm_update_loadings_scaled")
  free_idx <- setdiff(seq_len(setting$p), setting$fixed_loading_index)
  state$loading <- cpp_lgpjm_update_loadings_scaled(state$loading, state$y_star,
    state$omega, A_d, state$psi, as.integer(setting$loading_group), as.integer(free_idx),
    prior_var, prior_mean)
  state$loading[setting$fixed_loading_index] <- 1
  state
}

update_psi_continuous <- function(state, data, setting, a0 = 1, b0 = 1, A_d = NULL) {
  if (exists("cpp_update_psi", mode = "function")) {
    n_ord <- length(setting$ordinal_categories)
    state$psi <- cpp_update_psi(state$psi, state$y_star, state$omega, state$loading, A_d,
      as.integer(setting$loading_group), n_ord, a0, b0)
    state$psi[seq_len(n_ord)] <- 1
    return(state)
  }
  n <- dim(data$z_ord_obs)[3]
  Tn <- dim(data$z_ord_obs)[2]
  n_ord <- length(setting$ordinal_categories)
  for (j in (n_ord + 1):setting$p) {
    k <- setting$loading_group[j]
    rss <- 0
    cnt <- 0
    for (i in seq_len(n)) for (m in seq_len(Tn)) {
      mu <- A_d[j, i] + state$loading[j] * state$omega[k, m, i]
      rss <- rss + (state$y_star[j, m, i] - mu)^2
      cnt <- cnt + 1
    }
    state$psi[j] <- 1 / rgamma(1, shape = a0 + cnt / 2, rate = b0 + rss / 2)
  }
  state$psi[seq_len(n_ord)] <- 1
  state
}

update_latent_measurements <- function(state, data, setting, A_d, internal_method) {
  require_lgpjm_cpp("cpp_update_latent_measurements")
  n_ord <- length(setting$ordinal_categories)
  state$y_star <- cpp_update_latent_measurements(state$y_star, data$z_ord_obs_ref,
    data$x_cont_obs_ref, state$omega, state$loading, state$psi, A_d,
    state$thresholds, as.integer(setting$loading_group), n_ord, internal_method)
  state
}

update_thresholds_and_y <- function(state, data, setting, fit_cfg, A_d) {
  require_lgpjm_cpp("cpp_lgpjm_threshold_y_mh_mat")
  tm <- thresholds_to_matrix(state$thresholds)
  out <- cpp_lgpjm_threshold_y_mh_mat(tm$mat, tm$len, state$y_star, data$z_ord_obs_ref,
    state$omega, state$loading, A_d, as.integer(setting$loading_group),
    fit_cfg$threshold_proposal_var, state$accept_threshold, state$total_threshold)
  state$thresholds <- thresholds_from_matrix(out$thresholds_mat, tm$len)
  state$y_star <- out$y_star
  state$accept_threshold <- out$accept_threshold
  state$total_threshold <- out$total_threshold
  state
}

update_omega_mh <- function(state, data, setting, fit_cfg, A_d) {
  require_lgpjm_cpp("cpp_lgpjm_omega_mh")
  out <- cpp_lgpjm_omega_mh(state$omega, state$kesi_p, state$y_star, state$loading, state$psi,
    A_d, state$Qnew_inv, state$B_diag, state$alpha_ph_p, state$beta, state$lambda,
    data$obs_time, as.integer(data$delta), data$u_cov, state$cuts, setting$time_points,
    as.integer(setting$loading_group), as.integer(setting$N_int), setting$t_start,
    setting$sigmaf_omega, rep(setting$l_omega, setting$q), fit_cfg$gp_pred_nugget, fit_cfg$c_omega,
    state$accept_omega, state$total_omega)
  state$omega <- out$omega
  state$kesi_p <- out$kesi_p
  state$accept_omega <- out$accept_omega
  state$total_omega <- out$total_omega
  state
}

update_omega_linear <- function(state, data, setting, fit_cfg, A_d) {
  require_lgpjm_cpp("cpp_lgpjm_omega_linear_mh")
  out <- cpp_lgpjm_omega_linear_mh(state$omega, state$kesi_p, state$y_star, state$loading,
    state$psi, A_d, data$obs_time, as.integer(data$delta), data$u_cov, state$cuts,
    setting$time_points, state$alpha_ph_p, state$beta, state$lambda,
    as.integer(setting$loading_group), as.integer(setting$N_int), setting$t_start,
    fit_cfg$omega_linear_prior_var, fit_cfg$c_omega_linear, state$accept_omega, state$total_omega)
  state$omega <- out$omega
  state$kesi_p <- out$kesi_p
  state$accept_omega <- out$accept_omega
  state$total_omega <- out$total_omega
  state
}

update_alpha_ph_mh <- function(state, data, setting, fit_cfg) {
  require_lgpjm_cpp("cpp_lgpjm_alpha_ph_mh")
  out <- cpp_lgpjm_alpha_ph_mh(state$alpha_ph, state$alpha_ph_p, state$omega,
    state$kesi_p, state$beta, state$lambda, data$obs_time, as.integer(data$delta), data$u_cov, state$cuts,
    setting$time_points, as.integer(setting$N_int), setting$t_start, setting$sigma_alpha_ph,
    setting$l_alpha_ph, fit_cfg$gp_pred_nugget, setting$sigma_alpha_ph, setting$l_alpha_ph,
    fit_cfg$c_alpha, state$accept_alpha, state$total_alpha)
  state$alpha_ph <- out$alpha_ph
  state$alpha_ph_p <- out$alpha_ph_p
  state$accept_alpha <- out$accept_alpha
  state$total_alpha <- out$total_alpha
  state
}

update_beta_mh <- function(state, data, setting, fit_cfg) {
  require_lgpjm_cpp("cpp_lgpjm_beta_mh")
  out <- cpp_lgpjm_beta_mh(state$beta, state$alpha_ph_p, state$kesi_p, state$lambda,
    data$obs_time, as.integer(data$delta), data$u_cov, state$cuts, fit_cfg$beta_prior_mean,
    as.integer(setting$N_int), setting$t_start, fit_cfg$theta_prior_sd, fit_cfg$c_beta,
    state$accept_beta, state$total_beta)
  state$beta <- out$beta
  state$accept_beta <- out$accept_beta
  state$total_beta <- out$total_beta
  state
}

update_lambda_piecewise <- function(state, data, setting, fit_cfg) {
  require_lgpjm_cpp("cpp_lgpjm_lambda_gibbs")
  state$lambda <- cpp_lgpjm_lambda_gibbs(state$lambda, state$alpha_ph_p, state$kesi_p,
    state$beta, data$obs_time, as.integer(data$delta), data$u_cov, state$cuts,
    as.integer(setting$N_int), setting$t_start, fit_cfg$lambda_prior_shape,
    fit_cfg$lambda_prior_rate)
  state
}

compute_survival_terms <- function(state, data, setting) {
  require_lgpjm_cpp("cpp_lgpjm_survival_terms")
  cpp_lgpjm_survival_terms(state$alpha_ph_p, state$kesi_p, state$beta, state$lambda,
    data$obs_time, as.integer(data$delta), data$u_cov, state$cuts,
    as.integer(setting$N_int), setting$t_start)
}

compute_longitudinal_ll_by_subject <- function(state, data, setting, A_d = NULL) {
  if (exists("cpp_compute_longitudinal_ll", mode = "function")) {
    n_ord <- length(setting$ordinal_categories)
    return(cpp_compute_longitudinal_ll(state$omega, state$loading, state$psi, A_d,
      data$z_ord_obs_ref, data$x_cont_obs_ref, as.logical(data$obs_mask_ord_ref),
      as.logical(data$obs_mask_cont_ref), state$thresholds, as.integer(setting$loading_group), n_ord))
  }
  n <- dim(data$z_ord_obs_ref)[3]
  Tn <- dim(data$z_ord_obs_ref)[2]
  n_ord <- dim(data$z_ord_obs_ref)[1]
  n_cont <- dim(data$x_cont_obs_ref)[1]
  ll <- numeric(n)
  for (i in seq_len(n)) {
    s <- 0
    for (m in seq_len(Tn)) {
      for (j in seq_len(n_ord)) {
        if (!data$obs_mask_ord_ref[j, m, i]) next
        z <- data$z_ord_obs_ref[j, m, i]
        k <- setting$loading_group[j]
        mu <- A_d[j, i] + state$loading[j] * state$omega[k, m, i]
        bounds <- c(-Inf, state$thresholds[[j]], Inf)
        p_val <- pnorm(bounds[z + 1], mu, 1) - pnorm(bounds[z], mu, 1)
        s <- s + log(max(p_val, 1e-12))
      }
      for (j in seq_len(n_cont)) {
        if (!data$obs_mask_cont_ref[j, m, i]) next
        gj <- n_ord + j
        k <- setting$loading_group[gj]
        mu <- A_d[gj, i] + state$loading[gj] * state$omega[k, m, i]
        s <- s + dnorm(data$x_cont_obs_ref[j, m, i], mu, sqrt(state$psi[gj]), log = TRUE)
      }
    }
    ll[i] <- s
  }
  ll
}

run_single_chain <- function(data, setting, fit_cfg, chain_id, method = "internal",
                             seed = 123, print_progress = FALSE) {
  chain_seed <- seed + 1000 * chain_id
  state <- init_chain_state(data, setting, fit_cfg, chain_seed = chain_seed)
  state$cuts <- make_baseline_cuts(data$obs_time, data$delta, fit_cfg$G,
                                   fit_cfg$partition_method, setting$t_start, setting$t_end)

  estimate_gp_B <- isTRUE(fit_cfg$gp_estimate_sigma2 %||% TRUE) ||
    isTRUE(fit_cfg$gp_estimate_lengthscale %||% TRUE)
  free_idx <- setdiff(seq_len(setting$p), setting$fixed_loading_index)
  psi_idx <- (length(setting$ordinal_categories) + 1):setting$p
  gp_hyper_names <- if (fit_cfg$fit_model == "gp") {
    c(paste0("B", seq_len(setting$q)), paste0("l_B", seq_len(setting$q)))
  } else character()
  param_names <- c(
    paste0("beta", seq_along(state$beta)),
    paste0("alpha", seq_len(setting$q)),
    paste0("lambda", seq_along(state$lambda)),
    paste0("loading", free_idx),
    paste0("psi", psi_idx),
    gp_hyper_names
  )
  threshold_names <- threshold_param_names(state$thresholds)
  keep_n <- fit_cfg$n_iter - fit_cfg$burnin
  samples <- matrix(NA_real_, nrow = keep_n, ncol = length(param_names),
                    dimnames = list(NULL, param_names))
  store_chain_paths <- isTRUE(fit_cfg$save_chains)
  threshold_samples <- if (store_chain_paths && length(threshold_names) > 0) {
    matrix(NA_real_, nrow = keep_n, ncol = length(threshold_names),
           dimnames = list(NULL, threshold_names))
  } else NULL
  alpha_ph_samples <- if (store_chain_paths) {
    array(NA_real_, dim = c(keep_n, setting$q, setting$n_time),
          dimnames = list(NULL, paste0("alpha", seq_len(setting$q)),
                          paste0("t", seq_len(setting$n_time))))
  } else NULL
  loglik <- if (isTRUE(fit_cfg$compute_criteria)) {
    matrix(NA_real_, nrow = keep_n, ncol = setting$n_subject)
  } else NULL
  keep_id <- 0L

  for (iter in seq_len(fit_cfg$n_iter)) {
    if (print_progress && iter %% 1000 == 0) {
      cat(sprintf("[chain %d] iteration %d / %d\n", chain_id, iter, fit_cfg$n_iter))
      flush.console()
    }

    A_d <- compute_A_d(state, data)
    state <- update_latent_measurements(state, data, setting, A_d, internal_method = identical(method, "internal"))
    if (identical(fit_cfg$fit_model, "linear")) {
      state <- update_omega_linear(state, data, setting, fit_cfg, A_d)
    } else {
      state <- update_omega_mh(state, data, setting, fit_cfg, A_d)
    }

    state <- update_A(state, data, setting, fit_cfg)
    A_d <- compute_A_d(state, data)

    state <- update_loadings(state, data, setting, fit_cfg$loading_prior_var,
                             fit_cfg$loading_prior_mean %||% 0, A_d)

    if (estimate_gp_B && fit_cfg$fit_model == "gp") {
      state <- update_gp_hyperparams(state, data, setting, fit_cfg)
    }

    state <- update_psi_continuous(state, data, setting,
                                   fit_cfg$psi_prior_shape, fit_cfg$psi_prior_rate, A_d)

    state <- update_thresholds_and_y(state, data, setting, fit_cfg, A_d)

    state <- update_alpha_ph_mh(state, data, setting, fit_cfg)
    state <- update_beta_mh(state, data, setting, fit_cfg)
    state <- update_lambda_piecewise(state, data, setting, fit_cfg)

    if (iter > fit_cfg$burnin) {
      keep_id <- keep_id + 1L
      gp_hypers <- if (fit_cfg$fit_model == "gp") c(state$B_diag, state$l_omega) else numeric()
      samples[keep_id, ] <- c(state$beta, rowMeans(state$alpha_ph), state$lambda,
                              state$loading[free_idx], state$psi[psi_idx], gp_hypers)
      if (!is.null(threshold_samples)) {
        threshold_samples[keep_id, ] <- flatten_free_thresholds(state$thresholds)
      }
      if (!is.null(alpha_ph_samples)) {
        alpha_ph_samples[keep_id, , ] <- state$alpha_ph
      }
      if (isTRUE(fit_cfg$compute_criteria)) {
        terms_use <- compute_survival_terms(state, data, setting)
        ll_long <- compute_longitudinal_ll_by_subject(state, data, setting, A_d)
        loglik[keep_id, ] <- ll_long + terms_use$loglik
      }
    }
  }

  list(
    samples = samples,
    threshold_samples = threshold_samples,
    alpha_ph_samples = alpha_ph_samples,
    loglik = loglik,
    cuts = state$cuts,
    thresholds = state$thresholds,
    accept_omega_rate = state$accept_omega / max(1, state$total_omega),
    accept_alpha_rate = state$accept_alpha / max(1, state$total_alpha),
    accept_beta_rate = state$accept_beta / max(1, state$total_beta),
    accept_threshold_rate = state$accept_threshold / max(1, state$total_threshold),
    accept_l_rate = state$accept_l / max(1, state$total_l),
    accept_theta_rate = mean(c(state$accept_alpha / max(1, state$total_alpha),
                               state$accept_beta / max(1, state$total_beta)))
  )
}
