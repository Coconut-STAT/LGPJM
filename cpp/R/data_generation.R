sample_measurement_error <- function(mode, sd_target, is_continuous = TRUE) {
  mode <- tolower(mode %||% "normal")
  if (!is_continuous || mode == "normal") return(rnorm(1, mean = 0, sd = sd_target))
  if (mode == "t3_scaled") return(sd_target * rt(1, df = 3) / sqrt(3))
  if (mode == "shifted_gamma") {
    return(sd_target * (rgamma(1, shape = 5, rate = 4) - 5 / 4) / sqrt(5 / 16))
  }
  stop(sprintf("Unsupported error mode: %s", mode))
}


sample_measurement_errors_vec <- function(mode, sd_targets, is_continuous = rep(TRUE, length(sd_targets))) {
  n <- length(sd_targets)
  mode <- tolower(mode %||% "normal")
  out <- numeric(n)
  if (n == 0) return(out)
  if (length(is_continuous) != n) stop("is_continuous must have same length as sd_targets")
  is_continuous <- as.logical(is_continuous)
  if (any(!is_continuous)) {
    out[!is_continuous] <- rnorm(sum(!is_continuous), mean = 0, sd = sd_targets[!is_continuous])
  }
  if (mode == "normal") {
    out[is_continuous] <- rnorm(sum(is_continuous), mean = 0, sd = sd_targets[is_continuous])
    return(out)
  }
  # Keep the non-Gaussian shape while matching the target marker-specific
  # variance. Otherwise psi is correctly estimated around Var(epsilon), but
  # the analysis truth remains 0.3 and creates artificial positive bias.
  if (mode == "t3_scaled") {
    out[is_continuous] <- sd_targets[is_continuous] * rt(sum(is_continuous), df = 3) / sqrt(3)
    return(out)
  }
  if (mode == "shifted_gamma") {
    raw <- rgamma(sum(is_continuous), shape = 5, rate = 4) - 5 / 4
    out[is_continuous] <- sd_targets[is_continuous] * raw / sqrt(5 / 16)
    return(out)
  }
  stop(sprintf("Unsupported error mode: %s", mode))
}

sample_omega_trajectory <- function(time_points, mode, kernel_cfg) {
  n_time <- length(time_points)
  if (mode == "sinusoidal") {
    r1 <- rnorm(2); r2 <- rnorm(2)
    w1 <- 2 * sin(r1[1] + r1[2] * time_points)
    w2 <- 2 * cos(r2[1] + r2[2] * time_points)
    return(rbind(w1, w2))
  }
  if (mode == "linear") {
    r1 <- rnorm(2); r2 <- rnorm(2)
    w1 <- r1[1] + r1[2] * time_points
    w2 <- r2[1] + r2[2] * time_points
    return(rbind(w1, w2))
  }
  if (mode == "gp_kernel") {
    K <- se_kernel(time_points, time_points, sigma2 = kernel_cfg$sigma2,
                   length_scale = kernel_cfg$length_scale) + diag(kernel_cfg$nugget, n_time)
    w1 <- MASS::mvrnorm(1, rep(0, n_time), K)
    w2 <- MASS::mvrnorm(1, rep(0, n_time), K)
    return(rbind(w1, w2))
  }
  stop(sprintf("Unsupported omega mode: %s", mode))
}

generate_survival_times <- function(omega, u_cov, setting) {
  n_subject <- dim(omega)[3]
  time_dense <- seq(setting$t_start, setting$t_end, length.out = 800)
  alpha_dense <- setting$alpha_fun(time_dense)
  baseline_dense <- setting$baseline_fun(time_dense)
  t_event <- rep(Inf, n_subject)

  for (i in seq_len(n_subject)) {
    omega_dense <- matrix(NA, nrow = setting$q, ncol = length(time_dense))
    for (k in seq_len(setting$q)) {
      omega_dense[k, ] <- approx(x = setting$time_points, y = omega[k, , i],
                                  xout = time_dense, rule = 2)$y
    }
    eta <- as.numeric(u_cov[i, ] %*% setting$beta_true) + colSums(alpha_dense * omega_dense)
    haz <- baseline_dense * exp(eta)
    cum_haz <- cumsum(c(0, (haz[-1] + haz[-length(haz)]) * diff(time_dense) / 2))
    target <- -log(runif(1))
    hit <- which(cum_haz >= target)[1]
    if (!is.na(hit)) t_event[i] <- time_dense[hit]
  }

  target <- setting$censoring_rate
  censor_rate_at <- function(cmax) {
    event_before_end <- is.finite(t_event) & t_event <= setting$t_end
    probs <- rep(1, n_subject)
    probs[event_before_end] <- ifelse(t_event[event_before_end] >= cmax, 1, t_event[event_before_end] / cmax)
    mean(probs)
  }
  finite_event <- t_event[is.finite(t_event)]
  hi_base <- if (length(finite_event) > 0) max(setting$t_end, max(finite_event) * 3) else setting$t_end * 3
  lo <- 1e-3; hi <- hi_base
  for (iter in seq_len(60)) {
    mid <- 0.5 * (lo + hi)
    if (censor_rate_at(mid) > target) lo <- mid else hi <- mid
  }
  cmax <- 0.5 * (lo + hi)
  c_time <- runif(n_subject, min = 0, max = cmax)
  obs_time <- pmin(t_event, c_time, setting$t_end)
  delta <- as.integer(is.finite(t_event) & t_event <= c_time & t_event <= setting$t_end)
  list(t_true = t_event, c_time = c_time, obs_time = obs_time, delta = delta)
}

extend_omega_at_obs_time <- function(omega, obs_time, setting) {
  q <- dim(omega)[1]
  Tn <- dim(omega)[2]
  n <- dim(omega)[3]
  out <- array(NA_real_, dim = c(q, Tn + 1, n))
  out[, seq_len(Tn), ] <- omega
  for (i in seq_len(n)) {
    for (k in seq_len(q)) {
      out[k, Tn + 1, i] <- approx(
        x = setting$time_points, y = omega[k, , i],
        xout = obs_time[i], rule = 2
      )$y
    }
  }
  out
}

apply_missingness <- function(z_ord, x_cont, y_latent, d_cov, setting) {
  n_ord <- dim(z_ord)[1]; n_cont <- dim(x_cont)[1]
  n_time <- dim(z_ord)[2]; n_subject <- dim(z_ord)[3]
  p <- n_ord + n_cont
  mask <- array(FALSE, dim = c(p, n_time, n_subject))
  mechanism <- toupper(setting$missing_mechanism)
  base_p <- setting$missing_rate
  if (is.na(base_p) || base_p < 0) base_p <- 0
  if (base_p > 1) base_p <- 1

  for (i in seq_len(n_subject)) {
    d2 <- scale(d_cov[, 2])[i]
    for (m in seq_len(n_time)) {
      for (j in seq_len(p)) {
        p_miss <- if (mechanism == "MCAR") {
          base_p
        } else {
          if (base_p <= 0) 0 else if (base_p >= 1) 1 else {
            prev <- if (m == 1) 0 else y_latent[j, m - 1, i]
            lp <- qlogis(base_p) + 0.5 * d2 + 0.4 * prev
            min(0.9, max(1e-6, expit(lp)))
          }
        }
        if (runif(1) < p_miss) mask[j, m, i] <- TRUE
      }
    }
  }

  block_hits <- 0L
  for (i in seq_len(n_subject)) {
    for (m in seq_len(n_time)) {
      if (runif(1) < setting$block_missing_prob) {
        k <- sample(seq_len(setting$q), size = 1)
        idx <- which(setting$loading_group == k)
        mask[idx, m, i] <- TRUE
        block_hits <- block_hits + 1L
      }
    }
  }
  enforce_block <- isTRUE(setting$enforce_block_missing) || (setting$block_missing_prob > 0 || base_p > 0)
  if (enforce_block && block_hits == 0L) {
    i <- sample(seq_len(n_subject), size = 1); m <- sample(seq_len(n_time), size = 1)
    k <- sample(seq_len(setting$q), size = 1); idx <- which(setting$loading_group == k)
    mask[idx, m, i] <- TRUE; block_hits <- 1L
  }

  z_miss <- z_ord; x_miss <- x_cont
  for (j in seq_len(n_ord)) {
    idx <- which(mask[j, , ], arr.ind = TRUE)
    if (nrow(idx) > 0) for (r in seq_len(nrow(idx))) z_miss[j, idx[r, 1], idx[r, 2]] <- NA
  }
  for (j in seq_len(n_cont)) {
    global_j <- n_ord + j
    idx <- which(mask[global_j, , ], arr.ind = TRUE)
    if (nrow(idx) > 0) for (r in seq_len(nrow(idx))) x_miss[j, idx[r, 1], idx[r, 2]] <- NA
  }
  list(z_ord = z_miss, x_cont = x_miss, mask_all = mask, block_hits = block_hits)
}

generate_dataset <- function(setting, seed = 123) {
  set.seed(seed)
  n <- setting$n_subject; Tn <- setting$n_time; p <- setting$p
  n_ord <- length(setting$ordinal_categories); n_cont <- p - n_ord
  d_cov <- cbind(rep(1, n), rnorm(n), rt(n, df = 5))
  u_cov <- matrix(rnorm(n * setting$h), nrow = n, ncol = setting$h)

  omega <- array(NA, dim = c(setting$q, Tn, n))
  for (i in seq_len(n)) {
    omega[, , i] <- sample_omega_trajectory(setting$time_points, setting$omega_mode, setting$omega_kernel)
  }

  surv <- generate_survival_times(omega, u_cov, setting)
  omega_aug <- extend_omega_at_obs_time(omega, surv$obs_time, setting)
  T_aug <- dim(omega_aug)[2]

  y_latent <- array(NA, dim = c(p, T_aug, n))
  psi_sd <- sqrt(setting$psi_true)
  A_d_mat <- setting$A_true %*% t(d_cov)
  is_continuous_marker <- seq_len(p) > n_ord
  for (i in seq_len(n)) {
    a_d_i <- A_d_mat[, i]
    for (m in seq_len(T_aug)) {
      mu <- a_d_i + as.vector(setting$loading_true %*% omega_aug[, m, i])
      errs <- sample_measurement_errors_vec(setting$error_mode, psi_sd, is_continuous_marker)
      y_latent[, m, i] <- mu + errs
    }
  }

  thresholds <- vector("list", n_ord)
  z_ord <- array(NA_integer_, dim = c(n_ord, T_aug, n))
  for (j in seq_len(n_ord)) {
    K <- setting$ordinal_categories[j]
    probs <- seq(0, 1, length.out = K + 1)
    th <- as.numeric(quantile(y_latent[j, , ], probs = probs[2:K], names = FALSE))
    thresholds[[j]] <- th
    bounds <- c(-Inf, th, Inf)
    for (i in seq_len(n)) {
      for (m in seq_len(T_aug)) {
        z_ord[j, m, i] <- findInterval(y_latent[j, m, i], bounds, rightmost.closed = TRUE)
      }
    }
  }

  x_cont <- y_latent[(n_ord + 1):p, , , drop = FALSE]
  miss <- apply_missingness(z_ord, x_cont, y_latent, d_cov, setting)

  list(
    setting_id = setting$setting_id, seed = seed, time_points = setting$time_points,
    d_cov = d_cov, u_cov = u_cov, omega_true = omega_aug, y_latent_true = y_latent,
    z_ord_true = z_ord, x_cont_true = x_cont,
    z_ord_obs = miss$z_ord, x_cont_obs = miss$x_cont, thresholds = thresholds,
    obs_time = surv$obs_time, delta = surv$delta, censor_time = surv$c_time,
    event_time_true = surv$t_true, mask_all = miss$mask_all, missing_block_hits = miss$block_hits
  )
}
