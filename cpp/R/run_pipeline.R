## run_pipeline.R — Streamlined: one Simulation run → one .RData file

get_fit_config <- function(setting, cli) {
  n_cov <- ncol(setting$A_true)  # number of covariates
  # A prior mean: default 0 (Prior I), or from cli override
  A_prior_mean_val <- as_num(cli$A_prior_mean, 0)
  A_prior_mean <- matrix(A_prior_mean_val, nrow = setting$p, ncol = n_cov)
  # A prior variance: default I (Prior I), or 10^4*I (Prior II)
  A_prior_var <- as_num(cli$A_prior_var, 1)
  Sigma_A0_inv <- diag(1 / A_prior_var, n_cov)
  # Beta prior mean
  beta_prior_mean_val <- as_num(cli$beta_prior_mean, 0)
  beta_prior_mean <- rep(beta_prior_mean_val, setting$h)
  # Loading prior mean (for Prior II: Lambda_k0 = 2)
  loading_prior_mean_val <- as_num(cli$loading_prior_mean, 0)

  list(
    n_iter = as_int(cli$iterations, 10000),
    burnin = as_int(cli$burnin, 6000),
    n_chains = as_int(cli$chains, 1),
    n_datasets = as_int(cli$datasets, 500),
    G = as_int(cli$G, setting$baseline_partition_G),
    partition_method = as.character(cli$partition_method %||% setting$baseline_partition_method),
    rw_beta_sd = as_num(cli$rw_beta_sd, 0.08),
    rw_alpha_sd = as_num(cli$rw_alpha_sd, 0.06),
    c_omega = as_num(cli$c_omega, setting$c_omega %||% 0.5),
    c_omega_linear = as_num(cli$c_omega_linear, setting$c_omega_linear %||% 1),
    c_alpha = as_num(cli$c_alpha, as_num(cli$rw_alpha_sd, setting$c_alpha %||% 0.08)),
    c_beta = as_num(cli$c_beta, as_num(cli$rw_beta_sd, setting$c_beta %||% 0.08)),
    gp_pred_nugget = as_num(cli$gp_pred_nugget, setting$gp_pred_nugget %||% 0.1),
    threshold_proposal_var = as_num(cli$sigma_tau, (setting$sigma_tau %||% 0.1) / setting$n_time),
    theta_prior_sd = as_num(cli$theta_prior_sd, 1),
    lambda_prior_shape = as_num(cli$lambda_a1, 1),
    lambda_prior_rate = as_num(cli$lambda_a2, 0.01),
    loading_prior_var = as_num(cli$loading_prior_var, 1),
    loading_prior_mean = loading_prior_mean_val,
    psi_prior_shape = as_num(cli$psi_a0, 1),
    psi_prior_rate = as_num(cli$psi_b0, 1),
    omega_linear_prior_var = as_num(cli$omega_linear_prior_var, 100),
    fit_model = as.character(cli$fit_model %||% "gp"),
    methods = unique(trimws(strsplit(cli$methods %||% "internal", ",", fixed = TRUE)[[1]])),
    n_cores_requested = max(1L, as_int(cli$n_cores, as_int(Sys.getenv("N_CORES", unset = "500"), 500))),
    max_cores = max(1L, as_int(cli$max_cores, as_int(Sys.getenv("MAX_CORES", unset = "510"), 510))),
    force_cores = as_flag(cli$force_cores, FALSE),
    parallel_enabled = as_flag(cli$parallel, TRUE),
    save_chains = as_flag(cli$save_chains, FALSE),
    compute_criteria = as_flag(cli$compute_criteria, FALSE),
    A_prior_mean = A_prior_mean,
    Sigma_A0_inv = Sigma_A0_inv,
    beta_prior_mean = beta_prior_mean,
    gp_estimate_sigma2 = as_flag(cli$gp_estimate_sigma2, TRUE),
    gp_sigma2_a0 = as_num(cli$gp_sigma2_a0, 1),
    gp_sigma2_b0 = as_num(cli$gp_sigma2_b0, 1),
    gp_estimate_lengthscale = as_flag(cli$gp_estimate_lengthscale, TRUE),
    gp_l_prop_width = as_num(cli$gp_l_prop_width, setting$l_B_proposal_width %||% 0.2),
    gp_l_lower = as_num(cli$gp_l_lower, setting$l_B_lower %||% 0.02),
    gp_l_upper = as_num(cli$gp_l_upper, setting$l_B_upper %||% 5)
  )
}


prepare_data_for_method <- function(data, setting, method) {
  if (method == "imputation") {
    data <- impute_dataset(data, setting)
  }
  out <- data
  out$z_ord_obs_ref <- data$z_ord_obs
  out$x_cont_obs_ref <- data$x_cont_obs
  out$obs_mask_ord_orig <- !is.na(data$z_ord_obs)
  out$obs_mask_cont_orig <- !is.na(data$x_cont_obs)
  out$obs_mask_ord_ref <- !is.na(data$z_ord_obs)
  out$obs_mask_cont_ref <- !is.na(data$x_cont_obs)
  out
}

resolve_fit_model_by_method <- function(method, fit_cfg) {
  if (method == "linear_model") return("linear")
  fit_cfg$fit_model
}

get_truth_by_param <- function(setting, cuts) {
  free_idx <- setdiff(seq_len(setting$p), setting$fixed_loading_index)
  truth <- c(
    setNames(setting$beta_true, paste0("beta", seq_along(setting$beta_true))),
    setNames(rowMeans(setting$alpha_fun(setting$time_points)), paste0("alpha", seq_len(setting$q)))
  )
  lambda_truth <- numeric(length(cuts) - 1)
  for (g in seq_along(lambda_truth)) {
    seg <- seq(cuts[g], cuts[g + 1], length.out = 100)
    lambda_truth[g] <- mean(setting$baseline_fun(seg))
  }
  truth <- c(truth, setNames(lambda_truth, paste0("lambda", seq_along(lambda_truth))))
  loading_truth <- numeric(length(free_idx))
  names(loading_truth) <- paste0("loading", free_idx)
  for (m in seq_along(free_idx)) {
    j <- free_idx[m]; k <- setting$loading_group[j]
    loading_truth[m] <- setting$loading_true[j, k]
  }
  truth <- c(truth, loading_truth)
  psi_idx <- (length(setting$ordinal_categories) + 1):setting$p
  truth <- c(truth, setNames(setting$psi_true[psi_idx], paste0("psi", psi_idx)))
  truth
}

compute_group_rmse <- function(est_named, truth_named, group_pattern) {
  common <- intersect(names(est_named)[grepl(group_pattern, names(est_named))],
                       names(truth_named)[grepl(group_pattern, names(truth_named))])
  if (length(common) == 0) return(NA_real_)
  sqrt(mean((est_named[common] - truth_named[common])^2))
}

resolve_parallel_cores <- function(fit_cfg, n_tasks) {
  if (!fit_cfg$parallel_enabled || n_tasks <= 1) return(1L)
  if (isTRUE(fit_cfg$force_cores)) {
    return(as.integer(max(1L, min(fit_cfg$n_cores_requested, n_tasks))))
  }
  detect <- parallel::detectCores(logical = TRUE)
  if (is.na(detect) || detect < 1) detect <- 1L
  as.integer(max(1L, min(fit_cfg$n_cores_requested, fit_cfg$max_cores, max(1L, detect - 1L), n_tasks)))
}

# Assemble result for a single (rep, method) group
assemble_rep_result <- function(chains, setting) {
  sample_all <- combine_chain_matrices(chains, "samples")
  loglik_available <- all(vapply(chains, function(x) is.matrix(x$loglik), logical(1)))
  crit <- if (loglik_available) {
    loglik_all <- combine_chain_matrices(chains, "loglik")
    compute_waic_dic(loglik_all)
  } else {
    list(DIC = NA_real_, WAIC = NA_real_, pD = NA_real_, pWAIC = NA_real_, lppd = NA_real_)
  }
  summary_tbl <- summarize_parameters(sample_all)
  rhat <- compute_rhat(chains)
  cuts <- chains[[1]]$cuts
  truth <- get_truth_by_param(setting, cuts)
  est <- setNames(summary_tbl$mean, summary_tbl$parameter)
  q025 <- setNames(summary_tbl$q025, summary_tbl$parameter)
  q975 <- setNames(summary_tbl$q975, summary_tbl$parameter)
  common <- intersect(names(est), names(truth))
  param_eval <- data.frame(
    parameter = common, truth = unname(truth[common]),
    estimate = unname(est[common]), q025 = unname(q025[common]), q975 = unname(q975[common]),
    bias = unname(est[common] - truth[common]),
    sq_error = unname((est[common] - truth[common])^2),
    covered = as.integer(q025[common] <= truth[common] & q975[common] >= truth[common]),
    stringsAsFactors = FALSE
  )
  list(
    criteria = crit, cuts = cuts, rhat = rhat,
    measurement_rmse = compute_group_rmse(est, truth, "^(loading|psi)"),
    cox_rmse = compute_group_rmse(est, truth, "^beta"),
    accept_theta = mean(vapply(chains, function(x) x$accept_theta_rate, numeric(1))),
    accept_omega = mean(vapply(chains, function(x) x$accept_omega_rate %||% NA_real_, numeric(1)), na.rm = TRUE),
    accept_alpha = mean(vapply(chains, function(x) x$accept_alpha_rate %||% NA_real_, numeric(1)), na.rm = TRUE),
    accept_beta = mean(vapply(chains, function(x) x$accept_beta_rate %||% NA_real_, numeric(1)), na.rm = TRUE),
    accept_threshold = mean(vapply(chains, function(x) x$accept_threshold_rate %||% NA_real_, numeric(1)), na.rm = TRUE),
    accept_l = mean(vapply(chains, function(x) x$accept_l_rate %||% NA_real_, numeric(1)), na.rm = TRUE),
    param_eval = param_eval, summary = summary_tbl,
    truth = truth
  )
}

run_setting_pipeline <- function(setting, cli) {
  fit_cfg <- get_fit_config(setting, cli)
  seed <- if (!is.null(cli$seed) && nzchar(cli$seed)) as.integer(cli$seed) else 31L

  n_ds <- fit_cfg$n_datasets
  methods <- fit_cfg$methods
  n_chains <- fit_cfg$n_chains
  tasks <- list()
  for (r in seq_len(n_ds)) for (meth in methods) for (ch in seq_len(n_chains)) {
    tasks[[length(tasks) + 1]] <- list(rep_id = r, method = meth, chain = ch)
  }

  chain_worker <- function(task) {
    r <- task$rep_id; meth <- task$method; ch <- task$chain
    data_raw <- generate_dataset(setting, seed = seed + r)
    data <- prepare_data_for_method(data_raw, setting, meth)
    fit_model_use <- resolve_fit_model_by_method(meth, fit_cfg)
    fit_cfg_local <- fit_cfg; fit_cfg_local$fit_model <- fit_model_use
    chain <- run_single_chain(data, setting, fit_cfg_local, chain_id = ch,
                               method = meth, seed = seed + r,
                               print_progress = (r == 1 && ch == 1))
    list(rep_id = r, method = meth, chain_id = ch, chain = chain)
  }

  n_cores <- resolve_parallel_cores(fit_cfg, length(tasks))
  cat(sprintf("Running %d tasks (%d reps × %d methods × %d chains) on %d cores ...\n",
              length(tasks), n_ds, length(methods), n_chains, n_cores))

  if (n_cores > 1) {
    raw <- parallel::mclapply(tasks, chain_worker, mc.cores = n_cores, mc.preschedule = FALSE,
                               mc.silent = TRUE, mc.set.seed = TRUE)
  } else {
    raw <- lapply(tasks, chain_worker)
  }
  failed <- vapply(raw, function(x) inherits(x, "try-error") || is.null(x$chain), logical(1))
  if (any(failed)) cat(sprintf("WARNING: %d / %d workers failed. Skipping failed results.\n", sum(failed), length(raw)))
  raw <- raw[!failed]

  results <- list()
  for (cr in raw) {
    key <- paste(cr$rep_id, cr$method, sep = "_")
    if (is.null(results[[key]])) results[[key]] <- list(rep_id = cr$rep_id, method = cr$method, chains = list())
    results[[key]]$chains[[length(results[[key]]$chains) + 1]] <- cr$chain
  }

  all_results <- list()
  main_rows <- list()
  param_rows <- list()
  chain_objects <- list()
  for (key in names(results)) {
    r <- results[[key]]
    res <- tryCatch(assemble_rep_result(r$chains, setting), error = function(e) NULL)
    if (!is.null(res)) {
      res$rep_id <- r$rep_id; res$method <- r$method; res$setting <- setting$setting_id %||% "unknown"
      fit_model_label <- resolve_fit_model_by_method(res$method, fit_cfg)
      all_results[[length(all_results) + 1]] <- res
      main_rows[[length(main_rows) + 1]] <- data.frame(
        setting = res$setting,
        replicate = res$rep_id,
        method = res$method,
        fit_model = fit_model_label,
        DIC = res$criteria$DIC,
        WAIC = res$criteria$WAIC,
        measurement_rms = res$measurement_rmse,
        cox_rms = res$cox_rmse,
        theta_accept_rate = res$accept_theta,
        omega_accept_rate = res$accept_omega,
        alpha_accept_rate = res$accept_alpha,
        beta_accept_rate = res$accept_beta,
        threshold_accept_rate = res$accept_threshold,
        lengthscale_accept_rate = res$accept_l,
        stringsAsFactors = FALSE
      )
      if (nrow(res$param_eval) > 0) {
        prm <- res$param_eval
        prm$setting <- res$setting
        prm$replicate <- res$rep_id
        prm$method <- res$method
        prm$fit_model <- fit_model_label
        param_rows[[length(param_rows) + 1]] <- prm
      }
      if (isTRUE(fit_cfg$save_chains)) {
        chain_objects[[paste(res$rep_id, res$method, sep = "|")]] <- r$chains
      }
    }
  }
  list(
    setting = setting,
    fit_config = fit_cfg,
    main_results = if (length(main_rows)) do.call(rbind, main_rows) else data.frame(),
    parameter_results = if (length(param_rows)) do.call(rbind, param_rows) else data.frame(),
    chain_objects = if (isTRUE(fit_cfg$save_chains)) chain_objects else list(),
    rep_results = all_results
  )
}
