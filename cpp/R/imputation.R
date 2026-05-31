fill_locf_bidir <- function(x, default_value) {
  n <- length(x)
  if (n == 0L) return(x)
  if (all(is.na(x))) return(rep(default_value, n))
  out <- x
  last <- NA_real_
  for (i in seq_len(n)) {
    if (is.na(out[i])) {
      if (!is.na(last)) out[i] <- last
    } else {
      last <- out[i]
    }
  }
  next_val <- NA_real_
  for (i in n:1) {
    if (is.na(out[i])) {
      if (!is.na(next_val)) out[i] <- next_val
    } else {
      next_val <- out[i]
    }
  }
  out[is.na(out)] <- default_value
  out
}

impute_series_continuous <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  med <- stats::median(x, na.rm = TRUE)
  fill_locf_bidir(as.numeric(x), med)
}

impute_series_ordinal <- function(x, K) {
  if (all(is.na(x))) return(rep(ceiling(K / 2), length(x)))
  med <- stats::median(as.numeric(x), na.rm = TRUE)
  out <- fill_locf_bidir(as.numeric(x), med)
  out <- pmin(K, pmax(1, round(out)))
  as.integer(out)
}

impute_dataset <- function(data, setting) {
  z_imp <- data$z_ord_obs; x_imp <- data$x_cont_obs
  n_ord <- dim(z_imp)[1]; n_time <- dim(z_imp)[2]
  n_subject <- dim(z_imp)[3]; n_cont <- dim(x_imp)[1]
  for (j in seq_len(n_ord)) {
    K <- setting$ordinal_categories[j]
    for (i in seq_len(n_subject)) z_imp[j, , i] <- impute_series_ordinal(z_imp[j, , i], K)
  }
  for (j in seq_len(n_cont)) {
    for (i in seq_len(n_subject)) x_imp[j, , i] <- impute_series_continuous(x_imp[j, , i])
  }
  data$z_ord_obs <- z_imp; data$x_cont_obs <- x_imp; data$imputed <- TRUE
  data
}
