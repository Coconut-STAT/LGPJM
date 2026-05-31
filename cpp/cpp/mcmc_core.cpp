// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Squared exponential kernel
// [[Rcpp::export]]
arma::mat cpp_se_kernel(const arma::vec& x1, const arma::vec& x2,
                        double sigma2, double length_scale) {
  int n1 = x1.n_elem, n2 = x2.n_elem;
  arma::mat K(n1, n2);
  double denom = 2.0 * std::max(length_scale, 1e-12);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      double d = x1(i) - x2(j);
      K(i, j) = sigma2 * std::exp(-(d * d) / denom);
    }
  }
  return K;
}

// Truncated normal sampling (vectorized for one draw)
static double rtruncnorm_one(double mean, double sd, double lower, double upper) {
  double pl = R::pnorm(lower, mean, sd, 1, 0);
  double pu = R::pnorm(upper, mean, sd, 1, 0);
  if (pu - pl < 1e-15) {
    if (std::isfinite(lower) && std::isfinite(upper)) return 0.5 * (lower + upper);
    if (std::isfinite(lower)) return lower + 0.1 * sd;
    if (std::isfinite(upper)) return upper - 0.1 * sd;
    return mean;
  }
  double u = R::runif(pl, pu);
  double val = R::qnorm(u, mean, sd, 1, 0);
  if (!std::isfinite(val)) {
    if (std::isfinite(lower) && std::isfinite(upper)) val = 0.5 * (lower + upper);
    else if (std::isfinite(lower)) val = lower + 0.1 * sd;
    else if (std::isfinite(upper)) val = upper - 0.1 * sd;
    else val = mean;
  }
  return val;
}

static bool lgpjm_ord_bounds(double z_val, const arma::vec& tau,
                             double& lo, double& up) {
  if (!std::isfinite(z_val)) return false;
  int K = (int)tau.n_elem + 1;
  int z = (int)std::floor(z_val + 0.5);
  if (z < 1) z = 1;
  if (z > K) z = K;
  lo = (z <= 1) ? R_NegInf : tau(z - 2);
  up = (z >= K) ? R_PosInf : tau(z - 1);
  return lo < up;
}

static arma::mat lgpjm_inv_spd(const arma::mat& A);

// Update latent measurements (ordinal + continuous)
// [[Rcpp::export]]
arma::cube cpp_update_latent_measurements(
    arma::cube y_star,
    const arma::cube& z_ord_obs,
    const arma::cube& x_cont_obs,
    const arma::cube& omega,
    const arma::vec& loading,
    const arma::vec& psi,
    const arma::mat& A_d,
    const Rcpp::List& thresholds,
    const arma::ivec& loading_group,
    int n_ord,
    bool internal_method) {

  int p = y_star.n_rows;
  int Tn = y_star.n_cols;
  int n = y_star.n_slices;
  int n_cont = p - n_ord;

  for (int i = 0; i < n; i++) {
    for (int m = 0; m < Tn; m++) {
      // Ordinal markers
      for (int j = 0; j < n_ord; j++) {
        int k = loading_group(j) - 1; // 0-indexed
        double mu = A_d(j, i) + loading(j) * omega(k, m, i);
        double z_val = z_ord_obs(j, m, i);
        if (std::isnan(z_val)) {
          if (internal_method) {
            y_star(j, m, i) = R::rnorm(mu, 1.0);
          }
        } else {
          arma::vec th = thresholds[j];
          double lo = R_NegInf, up = R_PosInf;
          if (lgpjm_ord_bounds(z_val, th, lo, up)) {
            y_star(j, m, i) = rtruncnorm_one(mu, 1.0, lo, up);
          }
        }
      }
      // Continuous markers
      for (int jc = 0; jc < n_cont; jc++) {
        int gj = n_ord + jc;
        int k = loading_group(gj) - 1;
        double mu = A_d(gj, i) + loading(gj) * omega(k, m, i);
        double x_val = x_cont_obs(jc, m, i);
        if (std::isnan(x_val)) {
          if (internal_method) {
            y_star(gj, m, i) = R::rnorm(mu, std::sqrt(psi(gj)));
          }
        } else {
          y_star(gj, m, i) = x_val;
        }
      }
    }
  }
  return y_star;
}

// Update omega (GP model) for all subjects
// [[Rcpp::export]]
arma::cube cpp_update_omega_gp(
    arma::cube omega,
    const arma::cube& y_star,
    const arma::vec& loading,
    const arma::vec& psi,
    const arma::mat& A_d,
    const arma::mat& K_inv,
    const arma::ivec& loading_group,
    int q) {

  int Tn = omega.n_cols;
  int n = omega.n_slices;
  int p = y_star.n_rows;

  for (int i = 0; i < n; i++) {
    for (int k = 0; k < q; k++) {
      arma::vec w_diag(Tn, fill::zeros);
      arma::vec b_vec(Tn, fill::zeros);

      for (int m = 0; m < Tn; m++) {
        for (int j = 0; j < p; j++) {
          if (loading_group(j) - 1 != k) continue;
          double y_adj = y_star(j, m, i) - A_d(j, i);
          double lj = loading(j);
          double psi_inv = 1.0 / psi(j);
          w_diag(m) += lj * lj * psi_inv;
          b_vec(m) += lj * y_adj * psi_inv;
        }
      }

      arma::mat prec = K_inv + arma::diagmat(w_diag);
      arma::mat cov_post = arma::inv_sympd(prec);
      arma::vec mean_post = cov_post * b_vec;
      arma::vec draw = arma::mvnrnd(mean_post, cov_post);
      omega.subcube(k, 0, i, k, Tn - 1, i) = draw.t();
    }
  }
  return omega;
}

// Update omega (linear model) for all subjects
// [[Rcpp::export]]
arma::cube cpp_update_omega_linear(
    arma::cube omega,
    const arma::cube& y_star,
    const arma::vec& loading,
    const arma::vec& psi,
    const arma::mat& A_d,
    const arma::vec& time_points,
    const arma::ivec& loading_group,
    int q,
    double prior_var_linear) {

  int Tn = omega.n_cols;
  int n = omega.n_slices;
  int p = y_star.n_rows;

  for (int i = 0; i < n; i++) {
    for (int k = 0; k < q; k++) {
      arma::mat den = arma::eye(2, 2) / prior_var_linear;
      arma::vec num(2, fill::zeros);

      for (int m = 0; m < Tn; m++) {
        arma::vec x_t = {1.0, time_points(m)};
        for (int j = 0; j < p; j++) {
          if (loading_group(j) - 1 != k) continue;
          double psi_inv = 1.0 / psi(j);
          double y_adj = y_star(j, m, i) - A_d(j, i);
          arma::vec x_jm = loading(j) * x_t;
          den += x_jm * x_jm.t() * psi_inv;
          num += x_jm * y_adj * psi_inv;
        }
      }

      arma::mat cov_post = arma::inv_sympd(den);
      arma::vec mean_post = cov_post * num;
      arma::vec ab = arma::mvnrnd(mean_post, cov_post);
      for (int m = 0; m < Tn; m++) {
        omega(k, m, i) = ab(0) + ab(1) * time_points(m);
      }
    }
  }
  return omega;
}

// Update loadings
// [[Rcpp::export]]
arma::vec cpp_update_loadings(
    arma::vec loading,
    const arma::cube& y_star,
    const arma::cube& omega,
    const arma::mat& A_d,
    const arma::vec& psi,
    const arma::ivec& loading_group,
    const arma::ivec& free_idx,
    double prior_var) {

  int Tn = y_star.n_cols;
  int n = y_star.n_slices;

  for (int fi = 0; fi < (int)free_idx.n_elem; fi++) {
    int j = free_idx(fi) - 1; // 0-indexed
    int k = loading_group(j) - 1;
    double num = 0.0;
    double den = 1.0 / prior_var;
    double psi_inv = 1.0 / psi(j);

    for (int i = 0; i < n; i++) {
      for (int m = 0; m < Tn; m++) {
        double w = omega(k, m, i);
        double y_adj = y_star(j, m, i) - A_d(j, i);
        num += w * y_adj * psi_inv;
        den += w * w * psi_inv;
      }
    }
    double var_post = 1.0 / den;
    double mean_post = var_post * num;
    loading(j) = R::rnorm(mean_post, std::sqrt(var_post));
  }
  return loading;
}

// Update psi for continuous markers
// [[Rcpp::export]]
arma::vec cpp_update_psi(
    arma::vec psi,
    const arma::cube& y_star,
    const arma::cube& omega,
    const arma::vec& loading,
    const arma::mat& A_d,
    const arma::ivec& loading_group,
    int n_ord,
    double a0, double b0) {

  int p = y_star.n_rows;
  int Tn = y_star.n_cols;
  int n = y_star.n_slices;

  for (int j = n_ord; j < p; j++) {
    int k = loading_group(j) - 1;
    double rss = 0.0;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
      for (int m = 0; m < Tn; m++) {
        double mu = A_d(j, i) + loading(j) * omega(k, m, i);
        double resid = y_star(j, m, i) - mu;
        rss += resid * resid;
        cnt++;
      }
    }
    double shape = a0 + cnt / 2.0;
    double rate = b0 + rss / 2.0;
    psi(j) = 1.0 / R::rgamma(shape, 1.0 / rate);
  }
  return psi;
}

// Scaled-conjugate Gibbs update for A rows, matching the supplement/original code:
// A_j,-1 | psi_j ~ N(mean, psi_j * Sigma_A), with the intercept held fixed.
// [[Rcpp::export]]
arma::mat cpp_lgpjm_update_A_scaled(
    arma::mat A,
    const arma::cube& y_star,
    const arma::cube& omega,
    const arma::vec& loading,
    const arma::vec& psi,
    const arma::mat& d_cov,
    const arma::ivec& loading_group,
    const arma::mat& A_prior_mean,
    const arma::mat& Sigma_A0_inv,
    const arma::vec& fixed_intercepts) {

  int p = y_star.n_rows;
  int Tn = y_star.n_cols;
  int n = y_star.n_slices;
  int n_cov = d_cov.n_cols;
  if (n_cov < 2) {
    A.col(0) = fixed_intercepts;
    return A;
  }

  arma::mat D = d_cov.cols(1, n_cov - 1);
  arma::mat DtD = D.t() * D;
  arma::mat H0_inv = Sigma_A0_inv.submat(1, 1, n_cov - 1, n_cov - 1);

  for (int j = 0; j < p; ++j) {
    int k = loading_group(j) - 1;
    arma::vec resid_sum(n, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
      double s = 0.0;
      for (int m = 0; m < Tn; ++m) {
        s += y_star(j, m, i) - fixed_intercepts(j) * d_cov(i, 0) -
          loading(j) * omega(k, m, i);
      }
      resid_sum(i) = s;
    }

    arma::mat prec_base = H0_inv + Tn * DtD;
    arma::mat cov_base = lgpjm_inv_spd(prec_base);
    arma::vec prior_mean = A_prior_mean.row(j).cols(1, n_cov - 1).t();
    arma::vec rhs = H0_inv * prior_mean + D.t() * resid_sum;
    arma::vec mean = cov_base * rhs;
    arma::mat cov = std::max(psi(j), 1e-12) * cov_base;
    arma::vec draw = arma::mvnrnd(mean, cov);
    for (int c = 1; c < n_cov; ++c) A(j, c) = draw(c - 1);
  }

  A.col(0) = fixed_intercepts;
  return A;
}

// Scaled-conjugate Gibbs update for free factor loadings.
// [[Rcpp::export]]
arma::vec cpp_lgpjm_update_loadings_scaled(
    arma::vec loading,
    const arma::cube& y_star,
    const arma::cube& omega,
    const arma::mat& A_d,
    const arma::vec& psi,
    const arma::ivec& loading_group,
    const arma::ivec& free_idx,
    double prior_var,
    double prior_mean) {

  int Tn = y_star.n_cols;
  int n = y_star.n_slices;
  double v0 = std::max(prior_var, 1e-12);

  for (int fi = 0; fi < (int)free_idx.n_elem; ++fi) {
    int j = free_idx(fi) - 1;
    int k = loading_group(j) - 1;
    double ww = 0.0;
    double wy = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int m = 0; m < Tn; ++m) {
        double w = omega(k, m, i);
        double y_adj = y_star(j, m, i) - A_d(j, i);
        ww += w * w;
        wy += w * y_adj;
      }
    }
    double den = ww + 1.0 / v0;
    double mean = (wy + prior_mean / v0) / std::max(den, 1e-12);
    double var = std::max(psi(j), 1e-12) / std::max(den, 1e-12);
    loading(j) = R::rnorm(mean, std::sqrt(var));
  }
  return loading;
}

// Update continuous measurements; ordinal latent y is handled by the joint threshold+y MH.
// [[Rcpp::export]]
arma::cube cpp_lgpjm_update_continuous_measurements(
    arma::cube y_star,
    const arma::cube& x_cont_obs,
    const arma::cube& omega,
    const arma::vec& loading,
    const arma::vec& psi,
    const arma::mat& A_d,
    const arma::ivec& loading_group,
    int n_ord) {

  int n_cont = x_cont_obs.n_rows;
  int Tn = x_cont_obs.n_cols;
  int n = x_cont_obs.n_slices;
  for (int jc = 0; jc < n_cont; ++jc) {
    int gj = n_ord + jc;
    int k = loading_group(gj) - 1;
    double sd = std::sqrt(std::max(psi(gj), 1e-12));
    for (int i = 0; i < n; ++i) {
      for (int m = 0; m < Tn; ++m) {
        double x = x_cont_obs(jc, m, i);
        if (std::isnan(x)) {
          double mu = A_d(gj, i) + loading(gj) * omega(k, m, i);
          y_star(gj, m, i) = R::rnorm(mu, sd);
        } else {
          y_star(gj, m, i) = x;
        }
      }
    }
  }
  return y_star;
}

// Compute survival exposure and event terms
// [[Rcpp::export]]
Rcpp::List cpp_compute_survival_terms(
    const arma::vec& beta,
    const arma::vec& alpha,
    const arma::cube& omega,
    const arma::vec& obs_time,
    const arma::ivec& delta,
    const arma::mat& u_cov,
    const arma::vec& time_points,
    const arma::vec& cuts,
    int q) {

  int n = obs_time.n_elem;
  int G = cuts.n_elem - 1;
  int Tn = time_points.n_elem;
  arma::mat expo(n, G, fill::zeros);
  arma::vec eta_event(n, fill::zeros);
  arma::ivec event_g(n, fill::zeros);
  event_g.fill(-1);

  // Pre-compute beta*u for all subjects
  arma::vec beta_u = u_cov * beta;

  for (int i = 0; i < n; i++) {
    double t_i = obs_time(i);

    for (int g = 0; g < G; g++) {
      double lo = cuts(g);
      double up = std::min(cuts(g + 1), t_i);
      if (up <= lo) continue;

      // 10-point trapezoidal grid
      const int n_grid = 10;
      double dt = (up - lo) / (n_grid - 1);
      double sum = 0.0;
      for (int gi = 0; gi < n_grid; gi++) {
        double tt = lo + gi * dt;
        // Interpolate omega at tt
        double eta = beta_u(i);
        for (int k = 0; k < q; k++) {
          // Linear interpolation in time_points
          double wv;
          if (tt <= time_points(0)) {
            wv = omega(k, 0, i);
          } else if (tt >= time_points(Tn - 1)) {
            wv = omega(k, Tn - 1, i);
          } else {
            int idx = 0;
            for (int tm = 0; tm < Tn - 1; tm++) {
              if (time_points(tm + 1) >= tt) { idx = tm; break; }
            }
            double frac = (tt - time_points(idx)) / (time_points(idx + 1) - time_points(idx));
            wv = (1.0 - frac) * omega(k, idx, i) + frac * omega(k, idx + 1, i);
          }
          eta += alpha(k) * wv;
        }
        double val = std::exp(eta);
        if (gi == 0 || gi == n_grid - 1)
          sum += val;
        else
          sum += 2.0 * val;
      }
      expo(i, g) = sum * dt / 2.0;
    }

    if (delta(i) == 1) {
      int eg = 0;
      for (int g = 0; g < G; g++) {
        if (t_i >= cuts(g) && t_i <= cuts(g + 1)) { eg = g; break; }
      }
      if (eg >= G) eg = G - 1;
      event_g(i) = eg;

      double eta = beta_u(i);
      for (int k = 0; k < q; k++) {
        double wv;
        if (t_i <= time_points(0)) {
          wv = omega(k, 0, i);
        } else if (t_i >= time_points(Tn - 1)) {
          wv = omega(k, Tn - 1, i);
        } else {
          int idx = 0;
          for (int tm = 0; tm < Tn - 1; tm++) {
            if (time_points(tm + 1) >= t_i) { idx = tm; break; }
          }
          double frac = (t_i - time_points(idx)) / (time_points(idx + 1) - time_points(idx));
          wv = (1.0 - frac) * omega(k, idx, i) + frac * omega(k, idx + 1, i);
        }
        eta += alpha(k) * wv;
      }
      eta_event(i) = eta;
    }
  }

  // R convention: event_g is 1-indexed, -1 means no event
  arma::ivec event_g_r = event_g + 1; // shift to 1-based

  return List::create(
    Named("exposure") = expo,
    Named("event_g") = event_g_r,
    Named("eta_event") = eta_event
  );
}

// Compute longitudinal log-likelihood by subject
// [[Rcpp::export]]
arma::vec cpp_compute_longitudinal_ll(
    const arma::cube& omega,
    const arma::vec& loading,
    const arma::vec& psi,
    const arma::mat& A_d,
    const arma::cube& z_ord_obs_ref,
    const arma::cube& x_cont_obs_ref,
    const Rcpp::LogicalVector& obs_mask_ord_ref_flat,
    const Rcpp::LogicalVector& obs_mask_cont_ref_flat,
    const Rcpp::List& thresholds,
    const arma::ivec& loading_group,
    int n_ord) {

  int p = A_d.n_rows;
  int n_cont = p - n_ord;
  int Tn = z_ord_obs_ref.n_cols;
  int n_subj = z_ord_obs_ref.n_slices;
  arma::vec ll(n_subj, fill::zeros);

  for (int i = 0; i < n_subj; i++) {
    double s = 0.0;
    for (int m = 0; m < Tn; m++) {
      for (int j = 0; j < n_ord; j++) {
        int flat_idx = j + n_ord * (m + Tn * i);
        if (!obs_mask_ord_ref_flat[flat_idx]) continue;
        int k = loading_group(j) - 1;
        double mu = A_d(j, i) + loading(j) * omega(k, m, i);
        arma::vec th = thresholds[j];
        double lo_t = R_NegInf, up_t = R_PosInf;
        if (!lgpjm_ord_bounds(z_ord_obs_ref(j, m, i), th, lo_t, up_t)) continue;
        double prob = R::pnorm(up_t, mu, 1.0, 1, 0) - R::pnorm(lo_t, mu, 1.0, 1, 0);
        s += std::log(std::max(prob, 1e-12));
      }
      for (int jc = 0; jc < n_cont; jc++) {
        int flat_idx = jc + n_cont * (m + Tn * i);
        if (!obs_mask_cont_ref_flat[flat_idx]) continue;
        int gj = n_ord + jc;
        int k = loading_group(gj) - 1;
        double mu = A_d(gj, i) + loading(gj) * omega(k, m, i);
        double x = x_cont_obs_ref(jc, m, i);
        s += R::dnorm(x, mu, std::sqrt(psi(gj)), 1);
      }
    }
    ll(i) = s;
  }
  return ll;
}

// ---------------------------------------------------------------------------
// Reference-style LGPJM sampler helpers.
// These functions keep the current project's 3D array layout:
//   omega      [q, T + 1, n], last column is the subject-specific R_i value
//   y_star     [p, T, n], longitudinal latent/continuous values
//   kesi_p     [q, N_int + 1, n], omega interpolated to [a, R_i]
//   alpha_ph   [q, T], Cox time-varying coefficients on observation grid
//   alpha_ph_p [q, N_int + 1, n], alpha_ph interpolated to [a, R_i]
// ---------------------------------------------------------------------------

static arma::mat lgpjm_se_kernel(const arma::vec& x1, const arma::vec& x2,
                                 double sigmaf, double l) {
  int n1 = x1.n_elem, n2 = x2.n_elem;
  arma::mat K(n1, n2);
  double denom = 2.0 * std::max(l, 1e-12);
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      double d = x1(i) - x2(j);
      K(i, j) = sigmaf * std::exp(-(d * d) / denom);
    }
  }
  return K;
}

static arma::mat lgpjm_build_Qnew(const arma::vec& time_points,
                                  double sigmaf,
                                  double l,
                                  double nugget) {
  int Tn = time_points.n_elem;
  arma::mat Q = lgpjm_se_kernel(time_points, time_points, sigmaf, l);
  arma::mat Qnew(Tn + 1, Tn + 1, arma::fill::zeros);
  Qnew.submat(0, 0, Tn - 1, Tn - 1) = Q;
  Qnew(Tn, Tn) = 1.0;
  Qnew.diag() += nugget;
  return Qnew;
}

static double lgpjm_logdet_spd(const arma::mat& A) {
  double val = 0.0;
  double sign = 0.0;
  arma::log_det(val, sign, 0.5 * (A + A.t()));
  if (sign > 0.0 && std::isfinite(val)) return val;
  arma::vec eigval = arma::eig_sym(0.5 * (A + A.t()));
  eigval.transform([](double x) { return std::log(std::max(x, 1e-12)); });
  return arma::sum(eigval);
}

static arma::mat lgpjm_se_alpha_kernel(const arma::vec& x1, const arma::vec& x2,
                                       double sigmaf, double l,
                                       double cutoff = 0.25) {
  arma::mat K = lgpjm_se_kernel(x1, x2, sigmaf, l);
  for (arma::uword i = 0; i < K.n_rows; ++i) {
    for (arma::uword j = 0; j < K.n_cols; ++j) {
      if (std::fabs(x1(i) - x2(j)) > cutoff) K(i, j) = 0.0;
    }
  }
  return K;
}

static arma::mat lgpjm_inv_spd(const arma::mat& A) {
  arma::mat As = 0.5 * (A + A.t());
  arma::mat invA;
  if (arma::inv_sympd(invA, As)) return invA;
  arma::mat Aj = As;
  for (int k = 0; k < 6; ++k) {
    Aj.diag() += std::pow(10.0, -8 + k);
    if (arma::inv_sympd(invA, Aj)) return invA;
  }
  return arma::pinv(As);
}

static arma::mat lgpjm_chol_lower(const arma::mat& A) {
  arma::mat As = 0.5 * (A + A.t());
  arma::mat L;
  if (arma::chol(L, As, "lower")) return L;
  arma::mat Aj = As;
  for (int k = 0; k < 6; ++k) {
    Aj.diag() += std::pow(10.0, -8 + k);
    if (arma::chol(L, Aj, "lower")) return L;
  }
  arma::mat eigvec;
  arma::vec eigval;
  arma::eig_sym(eigval, eigvec, As);
  eigval.transform([](double x) { return std::sqrt(std::max(x, 1e-10)); });
  return eigvec * arma::diagmat(eigval);
}

static double lgpjm_safe_exp(double x) {
  if (x > 700.0) return std::exp(700.0);
  if (x < -700.0) return std::exp(-700.0);
  return std::exp(x);
}

static int lgpjm_segment_for_time(double tt, const arma::vec& cuts) {
  int G = cuts.n_elem - 1;
  if (tt <= cuts(0)) return 0;
  for (int g = 0; g < G; ++g) {
    if (tt <= cuts(g + 1) + 1e-12) return g;
  }
  return G - 1;
}

static arma::vec lgpjm_grid(double a, double end, int N_int) {
  arma::vec out(N_int + 1);
  if (N_int <= 0) {
    out.set_size(1);
    out(0) = end;
    return out;
  }
  for (int i = 0; i <= N_int; ++i) {
    out(i) = a + (end - a) * (double)i / (double)N_int;
  }
  return out;
}

static arma::mat lgpjm_predict_omega_subject(const arma::mat& omega_i,
                                             double obs_time_i,
                                             const arma::vec& time_points,
                                             int N_int,
                                             double a,
                                             double sigmaf,
                                             const arma::vec& l_vec,
                                             double nugget) {
  int q = omega_i.n_rows;
  int Tn = time_points.n_elem;
  int Tall = Tn + 1;
  arma::vec t_pred = lgpjm_grid(a, obs_time_i, N_int);
  arma::vec t_obs(Tall);
  arma::mat omega_obs(q, Tall, arma::fill::zeros);

  int happen = Tn;
  for (int m = 0; m < Tn; ++m) {
    if (obs_time_i <= time_points(m)) {
      happen = m;
      break;
    }
  }
  for (int m = 0; m < happen; ++m) {
    t_obs(m) = time_points(m);
    omega_obs.col(m) = omega_i.col(m);
  }
  t_obs(happen) = obs_time_i;
  omega_obs.col(happen) = omega_i.col(Tn);
  for (int m = happen + 1; m < Tall; ++m) {
    t_obs(m) = time_points(m - 1);
    omega_obs.col(m) = omega_i.col(m - 1);
  }
  arma::mat pred(q, N_int + 1, arma::fill::zeros);
  for (int k = 0; k < q; ++k) {
    double lk = l_vec.n_elem == 1 ? l_vec(0) : l_vec(k);
    arma::mat K = lgpjm_se_kernel(t_obs, t_obs, sigmaf, lk);
    K.diag() += nugget;
    arma::mat K_inv = lgpjm_inv_spd(K);
    arma::mat K_obs_pred = lgpjm_se_kernel(t_obs, t_pred, sigmaf, lk);
    pred.row(k) = omega_obs.row(k) * K_inv * K_obs_pred;
  }
  return pred;
}

static arma::mat lgpjm_predict_alpha_subject(const arma::mat& alpha_ph,
                                             double obs_time_i,
                                             const arma::vec& time_points,
                                             int N_int,
                                             double a,
                                             double sigmaf,
                                             double l,
                                             double nugget) {
  arma::vec t_pred = lgpjm_grid(a, obs_time_i, N_int);
  arma::mat K = lgpjm_se_kernel(time_points, time_points, sigmaf, l);
  K.diag() += nugget;
  arma::mat K_inv = lgpjm_inv_spd(K);
  arma::mat K_obs_pred = lgpjm_se_kernel(time_points, t_pred, sigmaf, l);
  arma::mat pred(alpha_ph.n_rows, N_int + 1, arma::fill::zeros);
  for (arma::uword k = 0; k < alpha_ph.n_rows; ++k) {
    pred.row(k) = alpha_ph.row(k) * K_inv * K_obs_pred;
  }
  return pred;
}

static void lgpjm_add_exposure_piece(double t0, double t1, double e0, double e1,
                                     const arma::vec& cuts, arma::vec& exposure) {
  if (t1 <= t0) return;
  double left = t0;
  int guard = 0;
  while (left < t1 - 1e-12 && guard++ < 100) {
    int g = lgpjm_segment_for_time(left + 1e-12, cuts);
    double right = std::min(t1, cuts(g + 1));
    if (right <= left + 1e-12) right = t1;
    double frac_l = (left - t0) / (t1 - t0);
    double frac_r = (right - t0) / (t1 - t0);
    double el = e0 + frac_l * (e1 - e0);
    double er = e0 + frac_r * (e1 - e0);
    exposure(g) += 0.5 * (el + er) * (right - left);
    left = right;
  }
}

static arma::vec lgpjm_exposure_one(const arma::mat& alpha_p_i,
                                    const arma::mat& kesi_p_i,
                                    const arma::vec& beta,
                                    const arma::rowvec& u_i,
                                    double obs_time_i,
                                    const arma::vec& cuts,
                                    int N_int,
                                    double a) {
  int G = cuts.n_elem - 1;
  arma::vec exposure(G, arma::fill::zeros);
  arma::vec t_int = lgpjm_grid(a, obs_time_i, N_int);
  double beta_u = arma::dot(beta, u_i.t());
  arma::vec e(N_int + 1);
  for (int r = 0; r <= N_int; ++r) {
    double eta = beta_u + arma::dot(alpha_p_i.col(r), kesi_p_i.col(r));
    e(r) = lgpjm_safe_exp(eta);
  }
  for (int r = 0; r < N_int; ++r) {
    lgpjm_add_exposure_piece(t_int(r), t_int(r + 1), e(r), e(r + 1), cuts, exposure);
  }
  return exposure;
}

static double lgpjm_survival_ll_one(const arma::mat& alpha_p_i,
                                    const arma::mat& kesi_p_i,
                                    const arma::vec& beta,
                                    const arma::vec& lambda,
                                    const arma::rowvec& u_i,
                                    double obs_time_i,
                                    int delta_i,
                                    const arma::vec& cuts,
                                    int N_int,
                                    double a) {
  arma::vec exposure = lgpjm_exposure_one(alpha_p_i, kesi_p_i, beta, u_i,
                                          obs_time_i, cuts, N_int, a);
  double ll = -arma::dot(lambda, exposure);
  if (delta_i == 1) {
    int g = lgpjm_segment_for_time(obs_time_i, cuts);
    double eta_event = arma::dot(beta, u_i.t()) +
      arma::dot(alpha_p_i.col(N_int), kesi_p_i.col(N_int));
    ll += std::log(std::max(lambda(g), 1e-300)) + eta_event;
  }
  return ll;
}

static double lgpjm_alpha_prior_ll(const arma::mat& alpha_ph,
                                   const arma::mat& K_alpha_inv) {
  double out = 0.0;
  for (arma::uword k = 0; k < alpha_ph.n_rows; ++k) {
    arma::rowvec arow = alpha_ph.row(k);
    out += -0.5 * arma::as_scalar(arow * K_alpha_inv * arow.t());
  }
  return out;
}

static double lgpjm_beta_prior_ll(const arma::vec& beta,
                                  const arma::vec& mean,
                                  double sd) {
  arma::vec d = beta - mean;
  return -0.5 * arma::dot(d, d) / std::max(sd * sd, 1e-12);
}

static double lgpjm_omega_prior_ll_subject(const arma::mat& omega_i,
                                           const arma::cube& Qnew_inv,
                                           const arma::vec& B_diag) {
  int q = omega_i.n_rows;
  double out = 0.0;
  for (int k = 0; k < q; ++k) {
    arma::rowvec w = omega_i.row(k);
    double inv_b = 1.0 / std::max(B_diag(k), 1e-12);
    out += -0.5 * inv_b * arma::as_scalar(w * Qnew_inv.slice(k) * w.t());
  }
  return out;
}

static double lgpjm_measure_ll_subject(const arma::mat& omega_i,
                                       const arma::cube& y_star,
                                       const arma::vec& loading,
                                       const arma::vec& psi,
                                       const arma::mat& A_d,
                                       const arma::ivec& loading_group,
                                       int subj) {
  int p = y_star.n_rows;
  int Tn = y_star.n_cols;
  double ll = 0.0;
  for (int j = 0; j < p; ++j) {
    int k = loading_group(j) - 1;
    double sd = std::sqrt(std::max(psi(j), 1e-12));
    for (int m = 0; m < Tn; ++m) {
      double mu = A_d(j, subj) + loading(j) * omega_i(k, m);
      ll += R::dnorm(y_star(j, m, subj), mu, sd, 1);
    }
  }
  return ll;
}

static double lgpjm_surv_ll_all(const arma::cube& alpha_ph_p,
                                const arma::cube& kesi_p,
                                const arma::vec& beta,
                                const arma::vec& lambda,
                                const arma::vec& obs_time,
                                const arma::ivec& delta,
                                const arma::mat& u_cov,
                                const arma::vec& cuts,
                                int N_int,
                                double a) {
  int n = obs_time.n_elem;
  double out = 0.0;
  for (int i = 0; i < n; ++i) {
    out += lgpjm_survival_ll_one(alpha_ph_p.slice(i), kesi_p.slice(i),
                                 beta, lambda, u_cov.row(i), obs_time(i),
                                 delta(i), cuts, N_int, a);
  }
  return out;
}

static double lgpjm_ord_marker_loglik(const arma::vec& tau,
                                      const arma::cube& z_ord_obs,
                                      const arma::cube& omega,
                                      const arma::vec& loading,
                                      const arma::mat& A_d,
                                      const arma::ivec& loading_group,
                                      int marker) {
  int Tn = z_ord_obs.n_cols;
  int n = z_ord_obs.n_slices;
  int k = loading_group(marker) - 1;
  double out = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int m = 0; m < Tn; ++m) {
      double z_val = z_ord_obs(marker, m, i);
      if (std::isnan(z_val)) continue;
      double lo = R_NegInf, up = R_PosInf;
      if (!lgpjm_ord_bounds(z_val, tau, lo, up)) continue;
      double mu = A_d(marker, i) + loading(marker) * omega(k, m, i);
      double pr = R::pnorm(up, mu, 1.0, 1, 0) - R::pnorm(lo, mu, 1.0, 1, 0);
      out += std::log(std::max(pr, 1e-14));
    }
  }
  return out;
}

static double lgpjm_log_interval_prob(double upper, double lower, double mean, double sd) {
  double pu = R::pnorm(upper, mean, sd, 1, 0);
  double pl = R::pnorm(lower, mean, sd, 1, 0);
  return std::log(std::max(pu - pl, 1e-14));
}

static void lgpjm_sample_ordinal_y_marker(arma::cube& y_star,
                                          const arma::vec& tau,
                                          const arma::cube& z_ord_obs,
                                          const arma::cube& omega,
                                          const arma::vec& loading,
                                          const arma::mat& A_d,
                                          const arma::ivec& loading_group,
                                          int marker) {
  int Tn = z_ord_obs.n_cols;
  int n = z_ord_obs.n_slices;
  int k = loading_group(marker) - 1;
  for (int i = 0; i < n; ++i) {
    for (int m = 0; m < Tn; ++m) {
      double z_val = z_ord_obs(marker, m, i);
      if (std::isnan(z_val)) continue;
      double lo = R_NegInf, up = R_PosInf;
      if (!lgpjm_ord_bounds(z_val, tau, lo, up)) continue;
      double mu = A_d(marker, i) + loading(marker) * omega(k, m, i);
      y_star(marker, m, i) = rtruncnorm_one(mu, 1.0, lo, up);
    }
  }
}

static arma::vec lgpjm_threshold_row_to_vec(const arma::mat& thresholds_mat,
                                            const arma::ivec& tau_len,
                                            int marker) {
  int max_len = thresholds_mat.n_cols;
  int s = tau_len(marker);
  if (s < 0) s = 0;
  if (s > max_len) s = max_len;
  arma::vec tau(s, arma::fill::zeros);
  for (int r = 0; r < s; ++r) tau(r) = thresholds_mat(marker, r);
  for (int r = 0; r < s; ++r) {
    if (!std::isfinite(tau(r))) tau(r) = (r == 0) ? 0.0 : tau(r - 1) + 0.25;
    if (r > 0 && tau(r) <= tau(r - 1)) tau(r) = tau(r - 1) + 0.25;
  }
  return tau;
}

// Matrix-based threshold + ordinal latent MH step.  This avoids repeatedly
// cloning and mutating an Rcpp::List inside heavily forked simulations.
// The first finite threshold for each ordinal marker is fixed for identifiability.
// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_threshold_y_mh_mat(arma::mat thresholds_mat,
                                        const arma::ivec& tau_len,
                                        arma::cube y_star,
                                        const arma::cube& z_ord_obs,
                                        const arma::cube& omega,
                                        const arma::vec& loading,
                                        const arma::mat& A_d,
                                        const arma::ivec& loading_group,
                                        double proposal_var,
                                        int accept_threshold,
                                        int total_threshold) {
  int n_ord = z_ord_obs.n_rows;
  double sd = std::sqrt(std::max(proposal_var, 1e-12));

  for (int j = 0; j < n_ord; ++j) {
    arma::vec tau = lgpjm_threshold_row_to_vec(thresholds_mat, tau_len, j);
    int s = tau.n_elem;
    for (int r = 0; r < s; ++r) thresholds_mat(j, r) = tau(r);
    if (s <= 1) {
      lgpjm_sample_ordinal_y_marker(y_star, tau, z_ord_obs, omega, loading,
                                    A_d, loading_group, j);
      continue;
    }

    total_threshold++;
    arma::vec prop = tau;
    for (int r = 1; r < s; ++r) {
      double lower = prop(r - 1);
      double upper = (r + 1 < s) ? tau(r + 1) : R_PosInf;
      prop(r) = rtruncnorm_one(tau(r), sd, lower, upper);
    }

    arma::vec old_ext(s + 2);
    arma::vec new_ext(s + 2);
    old_ext(0) = R_NegInf; new_ext(0) = R_NegInf;
    for (int r = 0; r < s; ++r) {
      old_ext(r + 1) = tau(r);
      new_ext(r + 1) = prop(r);
    }
    old_ext(s + 1) = R_PosInf; new_ext(s + 1) = R_PosInf;

    double log_prop_ratio = 0.0;
    bool ok = true;
    for (int h = 0; h < s - 1; ++h) {
      double num = R::pnorm((old_ext(h + 2) - old_ext(h + 1)) / sd, 0.0, 1.0, 1, 0) -
        R::pnorm((new_ext(h) - old_ext(h + 1)) / sd, 0.0, 1.0, 1, 0);
      double den = R::pnorm((new_ext(h + 2) - new_ext(h + 1)) / sd, 0.0, 1.0, 1, 0) -
        R::pnorm((old_ext(h) - new_ext(h + 1)) / sd, 0.0, 1.0, 1, 0);
      if (num <= 0.0 || den <= 0.0 || !std::isfinite(num) || !std::isfinite(den)) {
        ok = false;
        break;
      }
      log_prop_ratio += std::log(num) - std::log(den);
    }
    if (!ok) continue;

    double ll_old = lgpjm_ord_marker_loglik(tau, z_ord_obs, omega, loading,
                                            A_d, loading_group, j);
    double ll_new = lgpjm_ord_marker_loglik(prop, z_ord_obs, omega, loading,
                                            A_d, loading_group, j);
    if (std::log(R::runif(0.0, 1.0)) < ll_new - ll_old + log_prop_ratio) {
      accept_threshold++;
      for (int r = 0; r < s; ++r) thresholds_mat(j, r) = prop(r);
      lgpjm_sample_ordinal_y_marker(y_star, prop, z_ord_obs, omega, loading,
                                    A_d, loading_group, j);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("thresholds_mat") = thresholds_mat,
    Rcpp::Named("y_star") = y_star,
    Rcpp::Named("accept_threshold") = accept_threshold,
    Rcpp::Named("total_threshold") = total_threshold
  );
}

// Joint threshold + latent continuous ordinal measurement MH step.
// The first finite threshold for each ordinal marker is fixed for identifiability.
// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_threshold_y_mh(Rcpp::List thresholds,
                                    arma::cube y_star,
                                    const arma::cube& z_ord_obs,
                                    const arma::cube& omega,
                                    const arma::vec& loading,
                                    const arma::mat& A_d,
                                    const arma::ivec& loading_group,
                                    double proposal_var,
                                    int accept_threshold,
                                    int total_threshold) {
  int n_ord = z_ord_obs.n_rows;
  Rcpp::List out_thresholds = Rcpp::clone(thresholds);
  double sd = std::sqrt(std::max(proposal_var, 1e-12));

  for (int j = 0; j < n_ord; ++j) {
    arma::vec tau = Rcpp::as<arma::vec>(out_thresholds[j]);
    int s = tau.n_elem; // number of finite thresholds
    if (s <= 1) {
      lgpjm_sample_ordinal_y_marker(y_star, tau, z_ord_obs, omega, loading,
                                    A_d, loading_group, j);
      continue;
    }

    total_threshold++;
    arma::vec prop = tau;
    for (int r = 1; r < s; ++r) {
      double lower = prop(r - 1);
      double upper = (r + 1 < s) ? tau(r + 1) : R_PosInf;
      prop(r) = rtruncnorm_one(tau(r), sd, lower, upper);
    }

    arma::vec old_ext(s + 2);
    arma::vec new_ext(s + 2);
    old_ext(0) = R_NegInf; new_ext(0) = R_NegInf;
    for (int r = 0; r < s; ++r) {
      old_ext(r + 1) = tau(r);
      new_ext(r + 1) = prop(r);
    }
    old_ext(s + 1) = R_PosInf; new_ext(s + 1) = R_PosInf;

    double log_prop_ratio = 0.0;
    bool ok = true;
    for (int h = 0; h < s - 1; ++h) {
      double num = R::pnorm((old_ext(h + 2) - old_ext(h + 1)) / sd, 0.0, 1.0, 1, 0) -
        R::pnorm((new_ext(h) - old_ext(h + 1)) / sd, 0.0, 1.0, 1, 0);
      double den = R::pnorm((new_ext(h + 2) - new_ext(h + 1)) / sd, 0.0, 1.0, 1, 0) -
        R::pnorm((old_ext(h) - new_ext(h + 1)) / sd, 0.0, 1.0, 1, 0);
      if (num <= 0.0 || den <= 0.0 || !std::isfinite(num) || !std::isfinite(den)) {
        ok = false;
        break;
      }
      log_prop_ratio += std::log(num) - std::log(den);
    }
    if (!ok) continue;

    double ll_old = lgpjm_ord_marker_loglik(tau, z_ord_obs, omega, loading,
                                            A_d, loading_group, j);
    double ll_new = lgpjm_ord_marker_loglik(prop, z_ord_obs, omega, loading,
                                            A_d, loading_group, j);
    if (std::log(R::runif(0.0, 1.0)) < ll_new - ll_old + log_prop_ratio) {
      accept_threshold++;
      out_thresholds[j] = prop;
      lgpjm_sample_ordinal_y_marker(y_star, prop, z_ord_obs, omega, loading,
                                    A_d, loading_group, j);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("thresholds") = out_thresholds,
    Rcpp::Named("y_star") = y_star,
    Rcpp::Named("accept_threshold") = accept_threshold,
    Rcpp::Named("total_threshold") = total_threshold
  );
}

static arma::vec lgpjm_linear_fit_coef(const arma::vec& y_obs,
                                       const arma::vec& time_points) {
  int Tn = time_points.n_elem;
  arma::mat XtX(2, 2, arma::fill::zeros);
  arma::vec Xty(2, arma::fill::zeros);
  for (int m = 0; m < Tn; ++m) {
    double t = time_points(m);
    XtX(0, 0) += 1.0;
    XtX(0, 1) += t;
    XtX(1, 0) += t;
    XtX(1, 1) += t * t;
    Xty(0) += y_obs(m);
    Xty(1) += t * y_obs(m);
  }
  XtX.diag() += 1e-10;
  return lgpjm_inv_spd(XtX) * Xty;
}

static arma::rowvec lgpjm_eval_linear_values(const arma::vec& coef,
                                             const arma::vec& t_grid) {
  arma::rowvec out(t_grid.n_elem, arma::fill::zeros);
  for (arma::uword r = 0; r < t_grid.n_elem; ++r) {
    out(r) = coef(0) + coef(1) * t_grid(r);
  }
  return out;
}

// [[Rcpp::export]]
arma::cube cpp_lgpjm_compute_kesi_p(const arma::cube& omega,
                                    const arma::vec& obs_time,
                                    const arma::vec& time_points,
                                    int N_int,
                                    double a,
                                    double sigmaf,
                                    const arma::vec& l_vec,
                                    double nugget) {
  int q = omega.n_rows;
  int n = omega.n_slices;
  int Tall = omega.n_cols;
  arma::cube out(q, N_int + 1, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    arma::mat omega_i(q, Tall);
    for (int k = 0; k < q; ++k) {
      for (int m = 0; m < Tall; ++m) omega_i(k, m) = omega(k, m, i);
    }
    out.slice(i) = lgpjm_predict_omega_subject(omega_i, obs_time(i), time_points,
                                               N_int, a, sigmaf, l_vec, nugget);
  }
  return out;
}

// [[Rcpp::export]]
arma::cube cpp_lgpjm_compute_alpha_ph_p(const arma::mat& alpha_ph,
                                        const arma::vec& obs_time,
                                        const arma::vec& time_points,
                                        int N_int,
                                        double a,
                                        double sigmaf,
                                        double l,
                                        double nugget) {
  int q = alpha_ph.n_rows;
  int n = obs_time.n_elem;
  arma::cube out(q, N_int + 1, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    out.slice(i) = lgpjm_predict_alpha_subject(alpha_ph, obs_time(i), time_points,
                                               N_int, a, sigmaf, l, nugget);
  }
  return out;
}

// [[Rcpp::export]]
arma::cube cpp_lgpjm_compute_kesi_p_linear(const arma::cube& omega,
                                           const arma::vec& obs_time,
                                           const arma::vec& time_points,
                                           int N_int,
                                           double a) {
  int q = omega.n_rows;
  int n = omega.n_slices;
  int Tobs = time_points.n_elem;
  int usable_cols = std::min((int)omega.n_cols, Tobs + 1);
  arma::cube out(q, N_int + 1, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    arma::vec t_pred = lgpjm_grid(a, obs_time(i), N_int);
    arma::vec t_use(usable_cols, arma::fill::zeros);
    for (int m = 0; m < usable_cols; ++m) {
      t_use(m) = (m < Tobs) ? time_points(m) : obs_time(i);
    }
    for (int k = 0; k < q; ++k) {
      arma::vec y_obs(usable_cols, arma::fill::zeros);
      for (int m = 0; m < usable_cols; ++m) y_obs(m) = omega(k, m, i);
      arma::vec coef = lgpjm_linear_fit_coef(y_obs, t_use);
      out.slice(i).row(k) = lgpjm_eval_linear_values(coef, t_pred);
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_survival_terms(const arma::cube& alpha_ph_p,
                                    const arma::cube& kesi_p,
                                    const arma::vec& beta,
                                    const arma::vec& lambda,
                                    const arma::vec& obs_time,
                                    const arma::ivec& delta,
                                    const arma::mat& u_cov,
                                    const arma::vec& cuts,
                                    int N_int,
                                    double a) {
  int n = obs_time.n_elem;
  int G = cuts.n_elem - 1;
  arma::mat exposure(n, G, arma::fill::zeros);
  arma::vec eta_event(n, arma::fill::zeros);
  arma::ivec event_g(n);
  event_g.fill(NA_INTEGER);
  arma::vec loglik(n, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    exposure.row(i) = lgpjm_exposure_one(alpha_ph_p.slice(i), kesi_p.slice(i),
                                         beta, u_cov.row(i), obs_time(i),
                                         cuts, N_int, a).t();
    loglik(i) = -arma::dot(lambda, exposure.row(i).t());
    if (delta(i) == 1) {
      int g = lgpjm_segment_for_time(obs_time(i), cuts);
      event_g(i) = g + 1;
      eta_event(i) = arma::dot(beta, u_cov.row(i).t()) +
        arma::dot(alpha_ph_p.slice(i).col(N_int), kesi_p.slice(i).col(N_int));
      loglik(i) += std::log(std::max(lambda(g), 1e-300)) + eta_event(i);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("exposure") = exposure,
    Rcpp::Named("event_g") = event_g,
    Rcpp::Named("eta_event") = eta_event,
    Rcpp::Named("loglik") = loglik
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_omega_mh(arma::cube omega,
                              arma::cube kesi_p,
                              const arma::cube& y_star,
                              const arma::vec& loading,
                              const arma::vec& psi,
                              const arma::mat& A_d,
                              const arma::cube& Qnew_inv,
                              const arma::vec& B_diag,
                              const arma::cube& alpha_ph_p,
                              const arma::vec& beta,
                              const arma::vec& lambda,
                              const arma::vec& obs_time,
                              const arma::ivec& delta,
                              const arma::mat& u_cov,
                              const arma::vec& cuts,
                              const arma::vec& time_points,
                              const arma::ivec& loading_group,
                              int N_int,
                              double a,
                              double sigmaf,
                              const arma::vec& l_vec,
                              double nugget,
                              double c_omega,
                              int accept_omega,
                              int total_omega) {
  int q = omega.n_rows;
  int Tall = omega.n_cols;
  int Tn = y_star.n_cols;
  int n = omega.n_slices;
  int p = y_star.n_rows;

  std::vector<arma::mat> chol_props(q);
  for (int k = 0; k < q; ++k) {
    arma::vec w_diag(Tall, arma::fill::zeros);
    for (int m = 0; m < Tn; ++m) {
      for (int j = 0; j < p; ++j) {
        if (loading_group(j) - 1 != k) continue;
        w_diag(m) += loading(j) * loading(j) / std::max(psi(j), 1e-12);
      }
    }
    arma::mat prec = Qnew_inv.slice(k) / std::max(B_diag(k), 1e-12) + arma::diagmat(w_diag);
    arma::mat cov = lgpjm_inv_spd(prec);
    chol_props[k] = lgpjm_chol_lower(cov);
  }

  for (int i = 0; i < n; ++i) {
    arma::mat old_i(q, Tall);
    arma::mat prop_i(q, Tall);
    for (int k = 0; k < q; ++k) {
      arma::vec old_w(Tall);
      for (int m = 0; m < Tall; ++m) old_w(m) = omega(k, m, i);
      arma::vec prop_w = old_w + c_omega * chol_props[k] * arma::randn<arma::vec>(Tall);
      old_i.row(k) = old_w.t();
      prop_i.row(k) = prop_w.t();
    }

    arma::mat kesi_prop = lgpjm_predict_omega_subject(prop_i, obs_time(i), time_points,
                                                      N_int, a, sigmaf, l_vec, nugget);
    double lp_old = lgpjm_omega_prior_ll_subject(old_i, Qnew_inv, B_diag) +
      lgpjm_measure_ll_subject(old_i, y_star, loading, psi, A_d, loading_group, i) +
      lgpjm_survival_ll_one(alpha_ph_p.slice(i), kesi_p.slice(i), beta, lambda,
                            u_cov.row(i), obs_time(i), delta(i), cuts, N_int, a);
    double lp_new = lgpjm_omega_prior_ll_subject(prop_i, Qnew_inv, B_diag) +
      lgpjm_measure_ll_subject(prop_i, y_star, loading, psi, A_d, loading_group, i) +
      lgpjm_survival_ll_one(alpha_ph_p.slice(i), kesi_prop, beta, lambda,
                            u_cov.row(i), obs_time(i), delta(i), cuts, N_int, a);

    total_omega++;
    if (std::log(R::runif(0.0, 1.0)) < lp_new - lp_old) {
      accept_omega++;
      for (int k = 0; k < q; ++k) {
        for (int m = 0; m < Tall; ++m) omega(k, m, i) = prop_i(k, m);
      }
      kesi_p.slice(i) = kesi_prop;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("omega") = omega,
    Rcpp::Named("kesi_p") = kesi_p,
    Rcpp::Named("accept_omega") = accept_omega,
    Rcpp::Named("total_omega") = total_omega
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_omega_linear_mh(arma::cube omega,
                                     arma::cube kesi_p,
                                     const arma::cube& y_star,
                                     const arma::vec& loading,
                                     const arma::vec& psi,
                                     const arma::mat& A_d,
                                     const arma::vec& obs_time,
                                     const arma::ivec& delta,
                                     const arma::mat& u_cov,
                                     const arma::vec& cuts,
                                     const arma::vec& time_points,
                                     const arma::cube& alpha_ph_p,
                                     const arma::vec& beta,
                                     const arma::vec& lambda,
                                     const arma::ivec& loading_group,
                                     int N_int,
                                     double a,
                                     double prior_var_linear,
                                     double c_omega,
                                     int accept_omega,
                                     int total_omega) {
  int q = omega.n_rows;
  int Tobs = time_points.n_elem;
  int Tall = omega.n_cols;
  int Ty = y_star.n_cols;
  int Tfit = std::min(Tall, Ty);
  int n = omega.n_slices;
  int p = y_star.n_rows;
  double prior_sd2 = std::max(prior_var_linear, 1e-12);
  double c = std::isfinite(c_omega) && c_omega > 0.0 ? c_omega : 1.0;

  for (int i = 0; i < n; ++i) {
    arma::vec t_fit(Tfit, arma::fill::zeros);
    for (int m = 0; m < Tfit; ++m) {
      t_fit(m) = (m < Tobs) ? time_points(m) : obs_time(i);
    }
    arma::vec t_all(Tall, arma::fill::zeros);
    for (int m = 0; m < Tall; ++m) {
      t_all(m) = (m < Tobs) ? time_points(m) : obs_time(i);
    }

    arma::mat omega_old(q, Tfit, arma::fill::zeros);
    for (int k = 0; k < q; ++k) {
      for (int m = 0; m < Tfit; ++m) omega_old(k, m) = omega(k, m, i);
    }
    double surv_old = lgpjm_survival_ll_one(alpha_ph_p.slice(i), kesi_p.slice(i),
                                            beta, lambda, u_cov.row(i), obs_time(i),
                                            delta(i), cuts, N_int, a);
    arma::mat prop_i(q, Tall, arma::fill::zeros);
    arma::mat kesi_prop(q, N_int + 1, arma::fill::zeros);
    double qf_old_sum = 0.0;
    double qf_new_sum = 0.0;
    arma::vec t_grid = lgpjm_grid(a, obs_time(i), N_int);

    for (int k = 0; k < q; ++k) {
      arma::mat prec = arma::eye(2, 2) / prior_sd2;
      arma::vec num(2, arma::fill::zeros);
      for (int j = 0; j < p; ++j) {
        if (loading_group(j) - 1 != k) continue;
        double load = loading(j);
        double psi_inv = 1.0 / std::max(psi(j), 1e-12);
        for (int m = 0; m < Tfit; ++m) {
          double t = t_fit(m);
          double y_adj = y_star(j, m, i) - A_d(j, i);
          double w = load * load * psi_inv;
          prec(0, 0) += w;
          prec(0, 1) += w * t;
          prec(1, 0) += w * t;
          prec(1, 1) += w * t * t;
          num(0) += load * y_adj * psi_inv;
          num(1) += load * y_adj * t * psi_inv;
        }
      }

      arma::mat cov = lgpjm_inv_spd(prec);
      arma::mat L = lgpjm_chol_lower(cov);
      arma::vec mean = cov * num;
      arma::vec y_old = omega_old.row(k).t();
      arma::vec coef_old = lgpjm_linear_fit_coef(y_old, t_fit);
      arma::vec coef_new = mean + c * L * arma::randn<arma::vec>(2);

      arma::vec d_old = coef_old - mean;
      arma::vec d_new = coef_new - mean;
      qf_old_sum += arma::as_scalar(d_old.t() * prec * d_old);
      qf_new_sum += arma::as_scalar(d_new.t() * prec * d_new);

      for (int m = 0; m < Tall; ++m) {
        prop_i(k, m) = coef_new(0) + coef_new(1) * t_all(m);
      }
      kesi_prop.row(k) = lgpjm_eval_linear_values(coef_new, t_grid);
    }

    double surv_new = lgpjm_survival_ll_one(alpha_ph_p.slice(i), kesi_prop,
                                            beta, lambda, u_cov.row(i), obs_time(i),
                                            delta(i), cuts, N_int, a);
    double log_acc = surv_new - surv_old +
      0.5 * (1.0 - 1.0 / (c * c)) * (qf_old_sum - qf_new_sum);

    total_omega++;
    if (std::log(R::runif(0.0, 1.0)) < log_acc) {
      accept_omega++;
      for (int k = 0; k < q; ++k) {
        for (int m = 0; m < Tall; ++m) omega(k, m, i) = prop_i(k, m);
      }
      kesi_p.slice(i) = kesi_prop;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("omega") = omega,
    Rcpp::Named("kesi_p") = kesi_p,
    Rcpp::Named("accept_omega") = accept_omega,
    Rcpp::Named("total_omega") = total_omega
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_alpha_ph_mh(arma::mat alpha_ph,
                                 arma::cube alpha_ph_p,
                                 const arma::cube& omega,
                                 const arma::cube& kesi_p,
                                 const arma::vec& beta,
                                 const arma::vec& lambda,
                                 const arma::vec& obs_time,
                                 const arma::ivec& delta,
                                 const arma::mat& u_cov,
                                 const arma::vec& cuts,
                                 const arma::vec& time_points,
                                 int N_int,
                                 double a,
                                 double sigmaf,
                                 double l,
                                 double nugget,
                                 double sigma_alpha,
                                 double l_alpha,
                                 double c_alpha,
                                 int accept_alpha,
                                 int total_alpha) {
  arma::mat K_prior = lgpjm_se_kernel(time_points, time_points, sigma_alpha, l_alpha);
  K_prior.diag() += 1e-6;
  arma::mat K_prior_inv = lgpjm_inv_spd(K_prior);

  int q = alpha_ph.n_rows;
  int Tn = alpha_ph.n_cols;
  int n = obs_time.n_elem;
  arma::vec base_integral(n, arma::fill::zeros);
  arma::mat alpha_zero(q, N_int + 1, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    base_integral(i) = arma::dot(lambda, lgpjm_exposure_one(alpha_zero, kesi_p.slice(i),
                                                           beta, u_cov.row(i), obs_time(i),
                                                           cuts, N_int, a));
  }

  double lp_current = lgpjm_alpha_prior_ll(alpha_ph, K_prior_inv) +
    lgpjm_surv_ll_all(alpha_ph_p, kesi_p, beta, lambda, obs_time, delta,
                      u_cov, cuts, N_int, a);

  for (int k = 0; k < q; ++k) {
    arma::mat prec = K_prior_inv;
    for (int m = 0; m < Tn; ++m) {
      double info = 0.0;
      for (int i = 0; i < n; ++i) {
        double w = omega(k, m, i);
        info += base_integral(i) * w * w;
      }
      prec(m, m) += info;
    }
    arma::mat cov = lgpjm_inv_spd(prec);
    arma::mat L = lgpjm_chol_lower(cov);
    arma::vec step = c_alpha * L * arma::randn<arma::vec>(Tn);

    arma::mat alpha_new = alpha_ph;
    alpha_new.row(k) = alpha_new.row(k) + step.t();
    arma::cube alpha_p_new = cpp_lgpjm_compute_alpha_ph_p(alpha_new, obs_time, time_points,
                                                          N_int, a, sigmaf, l, nugget);
    double lp_new = lgpjm_alpha_prior_ll(alpha_new, K_prior_inv) +
      lgpjm_surv_ll_all(alpha_p_new, kesi_p, beta, lambda, obs_time, delta,
                        u_cov, cuts, N_int, a);

    total_alpha++;
    if (std::log(R::runif(0.0, 1.0)) < lp_new - lp_current) {
      accept_alpha++;
      alpha_ph = alpha_new;
      alpha_ph_p = alpha_p_new;
      lp_current = lp_new;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha_ph") = alpha_ph,
    Rcpp::Named("alpha_ph_p") = alpha_ph_p,
    Rcpp::Named("accept_alpha") = accept_alpha,
    Rcpp::Named("total_alpha") = total_alpha
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_beta_mh(arma::vec beta,
                             const arma::cube& alpha_ph_p,
                             const arma::cube& kesi_p,
                             const arma::vec& lambda,
                             const arma::vec& obs_time,
                             const arma::ivec& delta,
                             const arma::mat& u_cov,
                             const arma::vec& cuts,
                             const arma::vec& beta_prior_mean,
                             int N_int,
                             double a,
                             double beta_prior_sd,
                             double c_beta,
                             int accept_beta,
                             int total_beta) {
  int n = obs_time.n_elem;
  int h = beta.n_elem;
  arma::mat prec = arma::eye(h, h) / std::max(beta_prior_sd * beta_prior_sd, 1e-12);
  arma::vec beta_zero(h, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    double integral_no_beta = arma::dot(lambda, lgpjm_exposure_one(alpha_ph_p.slice(i),
                                                                   kesi_p.slice(i),
                                                                   beta_zero, u_cov.row(i),
                                                                   obs_time(i), cuts,
                                                                   N_int, a));
    arma::vec u = u_cov.row(i).t();
    prec += integral_no_beta * (u * u.t());
  }
  arma::mat cov = lgpjm_inv_spd(prec);
  arma::mat L = lgpjm_chol_lower(cov);
  arma::vec beta_new = beta + c_beta * L * arma::randn<arma::vec>(h);
  double lp_old = lgpjm_beta_prior_ll(beta, beta_prior_mean, beta_prior_sd) +
    lgpjm_surv_ll_all(alpha_ph_p, kesi_p, beta, lambda, obs_time, delta,
                      u_cov, cuts, N_int, a);
  double lp_new = lgpjm_beta_prior_ll(beta_new, beta_prior_mean, beta_prior_sd) +
    lgpjm_surv_ll_all(alpha_ph_p, kesi_p, beta_new, lambda, obs_time, delta,
                      u_cov, cuts, N_int, a);

  total_beta++;
  if (std::log(R::runif(0.0, 1.0)) < lp_new - lp_old) {
    accept_beta++;
    beta = beta_new;
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("accept_beta") = accept_beta,
    Rcpp::Named("total_beta") = total_beta
  );
}

// [[Rcpp::export]]
arma::vec cpp_lgpjm_lambda_gibbs(arma::vec lambda,
                                 const arma::cube& alpha_ph_p,
                                 const arma::cube& kesi_p,
                                 const arma::vec& beta,
                                 const arma::vec& obs_time,
                                 const arma::ivec& delta,
                                 const arma::mat& u_cov,
                                 const arma::vec& cuts,
                                 int N_int,
                                 double a,
                                 double prior_shape,
                                 double prior_rate) {
  int n = obs_time.n_elem;
  int G = cuts.n_elem - 1;
  arma::vec exposure_sum(G, arma::fill::zeros);
  arma::vec event_count(G, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    exposure_sum += lgpjm_exposure_one(alpha_ph_p.slice(i), kesi_p.slice(i),
                                       beta, u_cov.row(i), obs_time(i),
                                       cuts, N_int, a);
    if (delta(i) == 1) {
      event_count(lgpjm_segment_for_time(obs_time(i), cuts)) += 1.0;
    }
  }
  for (int g = 0; g < G; ++g) {
    double shape = prior_shape + event_count(g);
    double rate = prior_rate + exposure_sum(g);
    lambda(g) = R::rgamma(shape, 1.0 / std::max(rate, 1e-12));
  }
  return lambda;
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_B_gibbs(const arma::cube& omega,
                             const arma::mat& Qnew_inv,
                             arma::vec theta_phi,
                             double a0,
                             double b0) {
  int q = omega.n_rows;
  int Tall = omega.n_cols;
  int n = omega.n_slices;
  arma::vec B_diag(q, arma::fill::zeros);
  arma::vec quad(q, arma::fill::zeros);

  for (int k = 0; k < q; ++k) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
      arma::rowvec w(Tall);
      for (int m = 0; m < Tall; ++m) w(m) = omega(k, m, i);
      s += arma::as_scalar(w * Qnew_inv * w.t());
    }
    quad(k) = s;
    double shape = a0 + n * Tall / 2.0 + 1.0;
    double rate = b0 + s / 2.0;
    theta_phi(k) = R::rgamma(shape, 1.0 / std::max(rate, 1e-12));
  }
  if (q > 0) {
    for (int k = 0; k < q; ++k) {
      B_diag(k) = 1.0 / std::max(theta_phi(k), 1e-12);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("theta_phi") = theta_phi,
    Rcpp::Named("B_diag") = B_diag,
    Rcpp::Named("quad") = quad
  );
}

// Joint update for latent GP variance scales B and length-scales l_k.
// B is sampled by Gibbs; l_k uses the multiplicative MH proposal in the supplement.
// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_gp_hyper_update(const arma::cube& omega,
                                     arma::cube kesi_p,
                                     arma::cube Qnew_inv,
                                     arma::vec B_diag,
                                     arma::vec theta_phi,
                                     arma::vec l_vec,
                                     const arma::cube& alpha_ph_p,
                                     const arma::vec& beta,
                                     const arma::vec& lambda,
                                     const arma::vec& obs_time,
                                     const arma::ivec& delta,
                                     const arma::mat& u_cov,
                                     const arma::vec& cuts,
                                     const arma::vec& time_points,
                                     int N_int,
                                     double a,
                                     double sigmaf,
                                     double nugget,
                                     double a0,
                                     double b0,
                                     double l_prop_width,
                                     double l_lower,
                                     double l_upper,
                                     int accept_l,
                                     int total_l) {
  int q = omega.n_rows;
  int Tall = omega.n_cols;
  int n = omega.n_slices;
  double width = std::min(std::max(l_prop_width, 1e-6), 0.95);

  // Variance/precision Gibbs step.
  for (int k = 0; k < q; ++k) {
    double quad = 0.0;
    for (int i = 0; i < n; ++i) {
      arma::rowvec w(Tall);
      for (int m = 0; m < Tall; ++m) w(m) = omega(k, m, i);
      quad += arma::as_scalar(w * Qnew_inv.slice(k) * w.t());
    }
    double shape = a0 + n * Tall / 2.0 + 1.0;
    double rate = b0 + quad / 2.0;
    theta_phi(k) = R::rgamma(shape, 1.0 / std::max(rate, 1e-12));
    B_diag(k) = 1.0 / std::max(theta_phi(k), 1e-12);
  }

  // Length-scale MH step.
  for (int k = 0; k < q; ++k) {
    total_l++;
    double old_l = l_vec(k);
    double mult = R::runif(1.0 - width, 1.0 + width);
    double prop_l = old_l * mult;
    if (!(prop_l > l_lower && prop_l < l_upper)) continue;
    if (old_l < prop_l * (1.0 - width) || old_l > prop_l * (1.0 + width)) continue;

    arma::mat Qold = lgpjm_build_Qnew(time_points, sigmaf, old_l, nugget);
    arma::mat Qprop = lgpjm_build_Qnew(time_points, sigmaf, prop_l, nugget);
    arma::mat Qprop_inv = lgpjm_inv_spd(Qprop);
    double inv_b = 1.0 / std::max(B_diag(k), 1e-12);
    double quad_old = 0.0;
    double quad_new = 0.0;
    for (int i = 0; i < n; ++i) {
      arma::rowvec w(Tall);
      for (int m = 0; m < Tall; ++m) w(m) = omega(k, m, i);
      quad_old += arma::as_scalar(w * Qnew_inv.slice(k) * w.t());
      quad_new += arma::as_scalar(w * Qprop_inv * w.t());
    }
    double lp_old = -0.5 * n * lgpjm_logdet_spd(Qold) - 0.5 * inv_b * quad_old;
    double lp_new = -0.5 * n * lgpjm_logdet_spd(Qprop) - 0.5 * inv_b * quad_new;

    double log_hastings = std::log(std::max(old_l, 1e-12)) -
      std::log(std::max(prop_l, 1e-12));
    double log_acc = (lp_new - lp_old) + log_hastings;
    if (std::log(R::runif(0.0, 1.0)) < log_acc) {
      accept_l++;
      l_vec(k) = prop_l;
      Qnew_inv.slice(k) = Qprop_inv;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("theta_phi") = theta_phi,
    Rcpp::Named("B_diag") = B_diag,
    Rcpp::Named("l_vec") = l_vec,
    Rcpp::Named("Qnew_inv") = Qnew_inv,
    Rcpp::Named("kesi_p") = kesi_p,
    Rcpp::Named("accept_l") = accept_l,
    Rcpp::Named("total_l") = total_l
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_lgpjm_update_thresholds_rw(Rcpp::List thresholds,
                                          const arma::cube& z_ord_obs,
                                          const arma::cube& omega,
                                          const arma::vec& loading,
                                          const arma::mat& A_d,
                                          const arma::ivec& loading_group,
                                          double proposal_var,
                                          int accept_threshold,
                                          int total_threshold) {
  int n_ord = z_ord_obs.n_rows;
  Rcpp::List out_thresholds = Rcpp::clone(thresholds);
  double sd = std::sqrt(std::max(proposal_var, 1e-12));

  for (int j = 0; j < n_ord; ++j) {
    arma::vec tau = Rcpp::as<arma::vec>(out_thresholds[j]);
    int n_tau = tau.n_elem;
    if (n_tau <= 1) continue;

    total_threshold++;
    arma::vec prop = tau;
    for (int r = 1; r < n_tau; ++r) prop(r) = tau(r) + R::rnorm(0.0, sd);

    bool ordered = true;
    for (int r = 1; r < n_tau; ++r) {
      if (!(prop(r) > prop(r - 1))) {
        ordered = false;
        break;
      }
    }
    if (!ordered) continue;

    double ll_old = lgpjm_ord_marker_loglik(tau, z_ord_obs, omega, loading,
                                            A_d, loading_group, j);
    double ll_new = lgpjm_ord_marker_loglik(prop, z_ord_obs, omega, loading,
                                            A_d, loading_group, j);
    if (std::log(R::runif(0.0, 1.0)) < ll_new - ll_old) {
      accept_threshold++;
      out_thresholds[j] = prop;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("thresholds") = out_thresholds,
    Rcpp::Named("accept_threshold") = accept_threshold,
    Rcpp::Named("total_threshold") = total_threshold
  );
}
