`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.integer(x))
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x)) return(default)
  parts <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  vals <- suppressWarnings(as.numeric(trimws(parts)))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(default)
  vals
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      i <- i + 1
      next
    }
    key <- sub("^--", "", token)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

expit <- function(x) 1 / (1 + exp(-x))

trapz <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2) return(0)
  sum((y[-1] + y[-length(y)]) * diff(x) / 2)
}

se_kernel <- function(x1, x2, sigma2 = 1, length_scale = 1) {
  n1 <- length(x1); n2 <- length(x2)
  K <- matrix(0, n1, n2)
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      d <- x1[i] - x2[j]
      K[i, j] <- sigma2 * exp(-d^2 / (2 * length_scale))
    }
  }
  K
}

rtruncnorm1 <- function(mean, sd, lower, upper) {
  pl <- pnorm(lower, mean = mean, sd = sd)
  pu <- pnorm(upper, mean = mean, sd = sd)
  u <- runif(1, min = pl, max = pu)
  qnorm(u, mean = mean, sd = sd)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

timestamp_tag <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
