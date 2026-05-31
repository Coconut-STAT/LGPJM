combine_chain_matrices <- function(chains, field = "samples") {
  mats <- lapply(chains, function(x) x[[field]])
  do.call(rbind, mats)
}

compute_rhat <- function(chains) {
  mats <- lapply(chains, function(x) x$samples)
  m <- length(mats)
  if (m < 2) return(setNames(rep(NA_real_, ncol(mats[[1]])), colnames(mats[[1]])))
  n <- min(vapply(mats, nrow, integer(1)))
  mats <- lapply(mats, function(x) x[seq_len(n), , drop = FALSE])
  p <- ncol(mats[[1]])
  out <- numeric(p); names(out) <- colnames(mats[[1]])
  for (j in seq_len(p)) {
    chain_vals <- lapply(mats, function(x) x[, j])
    means <- vapply(chain_vals, mean, numeric(1))
    vars <- vapply(chain_vals, var, numeric(1))
    W <- mean(vars); B <- n * var(means)
    var_hat <- ((n - 1) / n) * W + (1 / n) * B
    out[j] <- sqrt(var_hat / W)
  }
  out
}

# Pretty label mapping for diagnostic traceplots
get_diagnostic_label <- function(param_name) {
  label_map <- c(
    "loading1" = expression(lambda[11]),
    "loading2" = expression(lambda[21]),
    "loading3" = expression(lambda[31]),
    "loading4" = expression(lambda[41]),
    "loading5" = expression(lambda[51]),
    "loading8" = expression(lambda[82]),
    "loading9" = expression(lambda[92]),
    "psi4"     = expression(psi[epsilon*4]),
    "psi5"     = expression(psi[epsilon*5]),
    "psi6"     = expression(psi[epsilon*6]),
    "psi7"     = expression(psi[epsilon*7]),
    "psi8"     = expression(psi[epsilon*8]),
    "psi9"     = expression(psi[epsilon*9]),
    "beta1"    = expression(beta[1]),
    "beta2"    = expression(beta[2]),
    "beta3"    = expression(beta[3]),
    "alpha1"   = expression(alpha[1](t)),
    "alpha2"   = expression(alpha[2](t))
  )
  if (param_name %in% names(label_map)) return(label_map[[param_name]])
  param_name
}

# Diagnostic traceplot: shows only loading, psi, beta with nice formatting
save_trace_plots <- function(chains, file_path, params_to_show = NULL) {
  mats <- lapply(chains, function(x) x$samples)
  param_names <- colnames(mats[[1]])

  if (is.null(params_to_show)) {
    # Default diagnostics focus on the reported measurement and Cox parameters.
    params_to_show <- param_names[grepl("^(loading|psi|beta)", param_names)]
  }
  keep_idx <- match(params_to_show, param_names)
  keep_idx <- keep_idx[!is.na(keep_idx)]
  n_p <- length(keep_idx)
  if (n_p == 0) return(invisible(NULL))

  # Layout: 4 columns, enough rows
  nc <- 4; nr <- ceiling(n_p / nc)

  # Color palette: elegant muted colors for chains
  chain_cols <- c("#2E86AB", "#E05263", "#4CAF50", "#FF9800", "#9C27B0")

  grDevices::pdf(file_path, width = 16, height = 3.2 * nr)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); grDevices::dev.off() }, add = TRUE)
  par(mfrow = c(nr, nc), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2, 0))

  for (pi in seq_along(keep_idx)) {
    idx <- keep_idx[pi]
    pname <- param_names[idx]
    max_it <- max(vapply(mats, nrow, integer(1)))
    all_vals <- unlist(lapply(mats, function(x) x[, idx]))
    val_range <- range(all_vals, na.rm = TRUE)
    # Expand y-axis by 15% for visual breathing room
    pad <- diff(val_range) * 0.15
    ylims <- c(val_range[1] - pad, val_range[2] + pad)

    plot(NA, xlim = c(1, max_it), ylim = ylims,
         xlab = "Iteration", ylab = "",
         main = "", axes = TRUE,
         cex.axis = 0.9, cex.lab = 1.0)
    title(main = get_diagnostic_label(pname), cex.main = 1.4)

    for (ch in seq_along(mats)) {
      lines(mats[[ch]][, idx], col = grDevices::adjustcolor(chain_cols[ch], alpha.f = 0.7), lwd = 0.8)
    }
  }

  # Fill empty panels
  if (n_p < nr * nc) {
    for (ei in seq_len(nr * nc - n_p)) plot.new()
  }
}

compute_waic_dic <- function(loglik_matrix) {
  if (!is.matrix(loglik_matrix)) stop("loglik_matrix must be a matrix")
  if (nrow(loglik_matrix) < 2) {
    return(list(DIC = NA_real_, WAIC = NA_real_, pD = NA_real_, pWAIC = NA_real_, lppd = NA_real_))
  }
  dev <- -2 * rowSums(loglik_matrix)
  Dbar <- mean(dev); pD <- stats::var(dev) / 2
  DIC <- Dbar + pD
  max_ll <- apply(loglik_matrix, 2, max)
  lppd <- sum(max_ll + log(colMeans(exp(sweep(loglik_matrix, 2, max_ll, "-")))))
  p_waic <- sum(apply(loglik_matrix, 2, var))
  WAIC <- -2 * (lppd - p_waic)
  list(DIC = DIC, WAIC = WAIC, pD = pD, pWAIC = p_waic, lppd = lppd)
}

summarize_parameters <- function(sample_matrix) {
  data.frame(
    parameter = colnames(sample_matrix),
    mean = colMeans(sample_matrix),
    sd = apply(sample_matrix, 2, sd),
    q025 = apply(sample_matrix, 2, stats::quantile, probs = 0.025),
    q50 = apply(sample_matrix, 2, stats::quantile, probs = 0.5),
    q975 = apply(sample_matrix, 2, stats::quantile, probs = 0.975),
    row.names = NULL
  )
}
