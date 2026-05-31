## mcmc_helpers.R — Utility functions for MCMC sampling
## SE kernels, GP regression, survival integrals

# SE kernel (standard)
se_kernel_mat <- function(x1, x2, sigmaf, l) {
  n1 <- length(x1); n2 <- length(x2)
  K <- matrix(0, n1, n2)
  for (i in seq_len(n1)) for (j in seq_len(n2)) {
    d <- x1[i] - x2[j]
    K[i, j] <- sigmaf * exp(-d^2 / (2 * l))
  }
  K
}

# SE kernel with cutoff (for alpha_ph GP prior)
se_kernel_alpha <- function(x1, x2, sigmaf, l, cutoff = 0.25) {
  n1 <- length(x1); n2 <- length(x2)
  K <- matrix(0, n1, n2)
  for (i in seq_len(n1)) for (j in seq_len(n2)) {
    d <- x1[i] - x2[j]
    K[i, j] <- sigmaf * exp(-d^2 / (2 * l))
    if (abs(d) > cutoff) K[i, j] <- 0
  }
  K
}

# Build Qnew: [Q, 0; 0, 1] — GP kernel extended with independent survival time point
build_Qnew <- function(Q, N) {
  Qnew <- cbind(Q, matrix(0, N, 1))
  Q1 <- matrix(0, 1, N + 1)
  Q1[1, N + 1] <- 1
  Qnew <- rbind(Qnew, Q1)
  Qnew
}

# GP regression: predict omega on integration grid for one subject
# t_obs: N+1 observed time points (reordered to include V_i)
# omega_obs: [q, N+1] omega values at t_obs
# t_pred: integration grid points
# sigmaf, l: kernel parameters
# Returns: [q, length(t_pred)] predicted omega values
gp_predict_omega <- function(omega_obs, t_obs, t_pred, sigmaf, l, nugget = 0.1) {
  q <- nrow(omega_obs)
  n_pred <- length(t_pred)
  result <- matrix(0, q, n_pred)
  K_obs <- se_kernel_mat(t_obs, t_obs, sigmaf, l)
  K_pred_obs <- se_kernel_mat(t_pred, t_obs, sigmaf, l)
  K_inv <- tryCatch(solve(K_obs + nugget * diag(length(t_obs))),
                    error = function(e) MASS::ginv(K_obs + nugget * diag(length(t_obs))))
  W <- K_pred_obs %*% K_inv  # [n_pred, N+1]
  for (k in seq_len(q)) {
    result[k, ] <- as.numeric(W %*% omega_obs[k, ])
  }
  result
}

# GP regression: predict alpha_ph on integration grid for one subject
# alpha_ph: [q, N] alpha values at observation times
# t_alpha: N observation time points
# t_pred: integration grid points
# Returns: [q, length(t_pred)]
gp_predict_alpha <- function(alpha_ph, t_alpha, t_pred, sigmaf, l, nugget = 0.1) {
  q <- nrow(alpha_ph)
  n_pred <- length(t_pred)
  result <- matrix(0, q, n_pred)
  K_obs <- se_kernel_mat(t_alpha, t_alpha, sigmaf, l)
  K_pred_obs <- se_kernel_mat(t_pred, t_alpha, sigmaf, l)
  K_inv <- tryCatch(solve(K_obs + nugget * diag(length(t_alpha))),
                    error = function(e) MASS::ginv(K_obs + nugget * diag(length(t_alpha))))
  W <- K_pred_obs %*% K_inv
  for (k in seq_len(q)) {
    result[k, ] <- as.numeric(W %*% alpha_ph[k, ])
  }
  result
}

# Compute kesi_p for all subjects: GP regression of omega onto integration grid
# omega: [q, N+1, n] — omega at N obs times + survival time
# V_i: [n] — observed survival/censoring times
# t_obs: [N] — regular observation time points
# N_int: number of integration intervals
# Returns: list with kesi_p [q, N_int+1, n]
compute_all_kesi_p <- function(omega, V_i, t_obs, N_int, a, sigmaf, l) {
  q <- dim(omega)[1]; n <- dim(omega)[3]
  N <- length(t_obs)
  kesi_p <- array(0, dim = c(q, N_int + 1, n))

  for (i in seq_len(n)) {
    # Integration grid from a to V_i[i]
    t_int <- seq(a, V_i[i], length.out = N_int + 1)

    # Build t_observed: insert V_i at correct position
    happen <- which(V_i[i] <= t_obs)[1]
    if (is.na(happen)) happen <- N + 1  # V_i beyond all obs times

    t_observed <- numeric(N + 1)
    omega_observed <- matrix(0, q, N + 1)

    if (happen > 1) {
      t_observed[1:(happen - 1)] <- t_obs[1:(happen - 1)]
      for (j in 1:(happen - 1)) omega_observed[, j] <- omega[, j, i]
    }
    t_observed[happen] <- V_i[i]
    omega_observed[, happen] <- omega[, N + 1, i]  # survival time value
    if (happen <= N) {
      t_observed[(happen + 1):(N + 1)] <- t_obs[happen:N]
      for (j in (happen + 1):(N + 1)) omega_observed[, j] <- omega[, j - 1, i]
    }

    kesi_p[, , i] <- gp_predict_omega(omega_observed, t_observed, t_int, sigmaf, l)
  }
  kesi_p
}

# Compute alpha_ph_p for all subjects
compute_all_alpha_ph_p <- function(alpha_ph, V_i, t_obs, N_int, a, sigmaf, l) {
  q <- nrow(alpha_ph); n <- length(V_i)
  alpha_ph_p <- array(0, dim = c(q, N_int + 1, n))
  for (i in seq_len(n)) {
    t_int <- seq(a, V_i[i], length.out = N_int + 1)
    alpha_ph_p[, , i] <- gp_predict_alpha(alpha_ph, t_obs, t_int, sigmaf, l)
  }
  alpha_ph_p
}

# Linear-trajectory interpolation of omega onto the survival integration grid
compute_all_kesi_p_linear <- function(omega, V_i, t_obs, N_int, a) {
  q <- dim(omega)[1]; n <- dim(omega)[3]
  Tn <- min(length(t_obs), dim(omega)[2])
  t_use <- t_obs[seq_len(Tn)]
  X <- cbind(1, t_use)
  XtX <- crossprod(X) + 1e-8 * diag(2)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
  kesi_p <- array(0, dim = c(q, N_int + 1, n))

  for (i in seq_len(n)) {
    t_int <- seq(a, V_i[i], length.out = N_int + 1)
    for (k in seq_len(q)) {
      y <- as.numeric(omega[k, seq_len(Tn), i])
      coef <- XtX_inv %*% crossprod(X, y)
      kesi_p[k, , i] <- as.numeric(coef[1] + coef[2] * t_int)
    }
  }
  kesi_p
}

# Compute survival integral for one subject using piecewise baseline hazard
# Matches reference f_integraltime.h
compute_integral_time_one <- function(alpha_ph_p_i, kesi_p_i, beta_ph, u_i,
                                       V_i_val, s, lambdak, N_int, a) {
  # alpha_ph_p_i: [q, N_int+1], kesi_p_i: [q, N_int+1]
  t_int <- seq(a, V_i_val, length.out = N_int + 1)
  integralresult <- 0
  flag_stop <- N_int + 1  # 1-indexed
  for (ii in seq_len(N_int + 1)) {
    if (V_i_val <= t_int[ii]) { flag_stop <- ii; break }
  }
  k <- 1  # current baseline segment (1-indexed)
  G <- length(lambdak)

  if (flag_stop > 2) {
    for (j in 1:(flag_stop - 2)) {
      eta_j <- sum(beta_ph * u_i) + sum(alpha_ph_p_i[, j] * kesi_p_i[, j])
      eta_j1 <- sum(beta_ph * u_i) + sum(alpha_ph_p_i[, j + 1] * kesi_p_i[, j + 1])

      if (k < G && t_int[j + 1] > s[k + 1]) {
        # Crosses segment boundary
        w1 <- (2 * t_int[j + 1] - s[k + 1] - t_int[j]) / (t_int[j + 1] - t_int[j])
        w2 <- (s[k + 1] - t_int[j]) / (t_int[j + 1] - t_int[j])
        c1 <- (exp(eta_j) * w1 + exp(eta_j1) * w2) / 2
        w3 <- (t_int[j + 1] - s[k + 1]) / (t_int[j + 1] - t_int[j])
        w4 <- (t_int[j + 1] + s[k + 1] - 2 * t_int[j]) / (t_int[j + 1] - t_int[j])
        c2 <- (exp(eta_j) * w3 + exp(eta_j1) * w4) / 2
        integralresult <- integralresult + lambdak[k] * c1 * (s[k + 1] - t_int[j]) +
          lambdak[k + 1] * c2 * (t_int[j + 1] - s[k + 1])
        k <- k + 1
      } else {
        cc <- (exp(eta_j) + exp(eta_j1)) / 2
        integralresult <- integralresult + lambdak[k] * cc * (t_int[j + 1] - t_int[j])
      }
    }
  }

  # Last segment to V_i
  eta_last <- sum(beta_ph * u_i) + sum(alpha_ph_p_i[, flag_stop - 1] * kesi_p_i[, flag_stop - 1])
  eta_end <- sum(beta_ph * u_i) + sum(alpha_ph_p_i[, flag_stop] * kesi_p_i[, flag_stop])

  if (k < G && V_i_val > s[k + 1]) {
    denom <- V_i_val - t_int[flag_stop - 1]
    if (denom < 1e-15) denom <- 1e-15
    w1 <- (2 * V_i_val - s[k + 1] - t_int[flag_stop - 1]) / denom
    w2 <- (s[k + 1] - t_int[flag_stop - 1]) / denom
    c1 <- (exp(eta_last) * w1 + exp(eta_end) * w2) / 2
    w3 <- (V_i_val - s[k + 1]) / denom
    w4 <- (V_i_val + s[k + 1] - 2 * t_int[flag_stop - 1]) / denom
    c2 <- (exp(eta_last) * w3 + exp(eta_end) * w4) / 2
    integralresult <- integralresult + lambdak[k] * c1 * (s[k + 1] - t_int[flag_stop - 1]) +
      lambdak[k + 1] * c2 * (V_i_val - s[k + 1])
  } else {
    cc <- (exp(eta_last) + exp(eta_end)) / 2
    integralresult <- integralresult + lambdak[k] * cc * (V_i_val - t_int[flag_stop - 1])
  }
  integralresult
}

# Truncated normal sampling (vectorized)
rtnorm_vec <- function(upper, lower, mean_vec, var, n_ignore) {
  sd_val <- sqrt(var)
  b <- pnorm(upper, mean = mean_vec, sd = sd_val)
  a_val <- pnorm(lower, mean = mean_vec, sd = sd_val)
  u <- runif(length(a_val), min = a_val, max = b)
  y <- qnorm(u, mean = mean_vec, sd = sd_val)
  # Handle Inf/-Inf
  idx_inf <- which(y == Inf)
  if (length(idx_inf) > 0) for (ii in idx_inf) y[ii] <- runif(1, lower[ii], lower[ii] + 0.1 * sd_val)
  idx_ninf <- which(y == -Inf)
  if (length(idx_ninf) > 0) for (ii in idx_ninf) y[ii] <- runif(1, upper[ii] - 0.1 * sd_val, upper[ii])
  idx_nan <- which(is.nan(y) | is.na(y))
  if (length(idx_nan) > 0) y[idx_nan] <- mean_vec[idx_nan]
  y
}
