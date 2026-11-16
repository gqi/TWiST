#include <Rcpp.h>

using namespace Rcpp;

// Group coordinate descent using individual-level summaries, with warm start along the lambda path
// [[Rcpp::export]]
List flasso_compressed_rep_cd_fit_cpp(
    NumericMatrix G,
    NumericMatrix B,
    IntegerVector donor_index,
    NumericVector y,
    NumericVector lambda,
    NumericVector feature_off,
    NumericVector Q_array,
    NumericVector group_multiplier,
    double alpha,
    double eps,
    int max_iter,
    NumericVector offset_vec) {
  const int D = G.nrow();
  const int J = G.ncol();
  const int n = B.nrow();
  const int K = B.ncol();
  const int p = J * K;
  const int L = lambda.size();

  if (donor_index.size() != n || y.size() != n || offset_vec.size() != n ||
      feature_off.size() != p || Q_array.size() != p * K ||
      group_multiplier.size() != J || L == 0)
    stop("Incompatible dimensions in the fast-TWiST solver inputs.");

  std::vector<int> donor(n);
  for (int i = 0; i < n; ++i) {
    donor[i] = donor_index[i] - 1;
    if (donor[i] < 0 || donor[i] >= D) stop("donor_index out of range");
  }

  const double* g_ptr = G.begin();
  const double* b_ptr = B.begin();
  const double* q_ptr = Q_array.begin();

  // S: sums of the B-spline bases within each individual (D x K)
  // M: sums of their cross-products (D x K x K)
  std::vector<double> S(D * K, 0.0), M(D * K * K, 0.0);
  for (int i = 0; i < n; ++i) {
    const int d = donor[i];
    for (int k = 0; k < K; ++k) {
      const double bik = b_ptr[i + n * k];
      S[d + D * k] += bik;
      for (int l = 0; l < K; ++l) M[d + D * (k + K * l)] += bik * b_ptr[i + n * l];
    }
  }

  NumericMatrix beta(p, L);
  NumericVector intercept(L);
  IntegerVector iter(L);
  LogicalVector converged(L);

  double ysum = 0.0;
  double exposure = 0.0;
  for (int i = 0; i < n; ++i) {
    ysum += y[i];
    exposure += std::exp(offset_vec[i]);
  }
  double a0 = std::log(ysum / exposure);

  std::vector<double> a(p, 0.0), bcoef(p, 0.0);
  std::vector<double> b_new(K), delta(K), z(K), q_delta(K), gh(K), eta(n);
  std::vector<double> H(D * K, 0.0), F(D * K, 0.0);
  std::vector<int> active(J, 0);
  double feature_offset = 0.0;
  double rsum = 0.0;
  int total_iter = 0;

  // Linear predictor of each cell, returns the maximum
  auto compute_eta = [&]() {
    double max_eta = R_NegInf;
    for (int i = 0; i < n; ++i) {
      double e = offset_vec[i] + a0 + feature_offset;
      for (int k = 0; k < K; ++k) e += b_ptr[i + n * k] * F[donor[i] + D * k];
      eta[i] = e;
      max_eta = std::max(max_eta, e);
    }
    return max_eta;
  };

  // Working residuals r and H = sum of B * r within each individual
  auto compute_H = [&](double v) {
    std::fill(H.begin(), H.end(), 0.0);
    rsum = 0.0;
    for (int i = 0; i < n; ++i) {
      const double ri = (y[i] - std::exp(eta[i])) / v;
      rsum += ri;
      for (int k = 0; k < K; ++k) H[donor[i] + D * k] += b_ptr[i + n * k] * ri;
    }
  };

  // Update F, H and sum(r) after the coefficients of SNP g change by del
  auto update_group = [&](int g, const std::vector<double>& del) {
    std::fill(q_delta.begin(), q_delta.end(), 0.0);
    double off_delta = 0.0;
    for (int h = 0; h < K; ++h) {
      off_delta += feature_off[g * K + h] * del[h];
      for (int l = 0; l < K; ++l) q_delta[l] += q_ptr[l + K * (h + K * g)] * del[h];
    }

    feature_offset += off_delta;
    double rsum_shift = off_delta * n;
    for (int l = 0; l < K; ++l) {
      const double qd = q_delta[l];
      if (qd == 0.0) continue;
      for (int d = 0; d < D; ++d) {
        const double gd = g_ptr[d + D * g];
        F[d + D * l] += gd * qd;
        rsum_shift += gd * S[d + D * l] * qd;
      }
    }
    rsum -= rsum_shift;

    if (off_delta != 0.0) {
      for (int k = 0; k < K; ++k)
        for (int d = 0; d < D; ++d) H[d + D * k] -= off_delta * S[d + D * k];
    }
    for (int l = 0; l < K; ++l) {
      const double qd = q_delta[l];
      if (qd == 0.0) continue;
      for (int k = 0; k < K; ++k) {
        for (int d = 0; d < D; ++d) {
          const double gd = g_ptr[d + D * g];
          H[d + D * k] -= gd * M[d + D * (k + K * l)] * qd;
        }
      }
    }
  };

  // Group update for SNP g, returns whether the group is nonzero
  auto fit_group = [&](int g, int ell, double v, double& max_change) {
    for (int l = 0; l < K; ++l) {
      double s = 0.0;
      for (int d = 0; d < D; ++d) s += g_ptr[d + D * g] * H[d + D * l];
      gh[l] = s;
    }

    for (int h = 0; h < K; ++h) {
      double cp = feature_off[g * K + h] * rsum;
      for (int l = 0; l < K; ++l) cp += q_ptr[l + K * (h + K * g)] * gh[l];
      z[h] = cp / n + a[g * K + h];
    }

    double z_ss = 0.0;
    for (int h = 0; h < K; ++h) z_ss += z[h] * z[h];
    const double z_norm = std::sqrt(z_ss);
    std::fill(b_new.begin(), b_new.end(), 0.0);
    const double mult = group_multiplier[g];
    const double lam1 = lambda[ell] * mult * alpha;
    const double lam2 = lambda[ell] * mult * (1.0 - alpha);
    if (z_norm > 0.0) {
      const double len = std::max(v * z_norm - lam1, 0.0) / (v * (1.0 + lam2));
      for (int h = 0; h < K; ++h) b_new[h] = len * z[h] / z_norm;
    }

    bool changed = false;
    for (int h = 0; h < K; ++h) {
      const int col = g * K + h;
      delta[h] = b_new[h] - a[col];
      if (delta[h] != 0.0) {
        changed = true;
        max_change = std::max(max_change, std::abs(delta[h]));
      }
      bcoef[col] = b_new[h];
    }
    if (changed) update_group(g, delta);
    return bcoef[g * K] != 0.0;
  };

  for (int ell = 0; ell < L; ++ell) {
    while (total_iter < max_iter) {
      bool active_converged = false;
      double v = 1.0;

      // Iterate over the active set until convergence
      while (total_iter < max_iter) {
        iter[ell]++;
        total_iter++;

        const double max_eta = compute_eta();
        v = std::exp(max_eta);
        compute_H(v);

        // Intercept
        const double shift0 = rsum / n;
        double max_change = std::abs(shift0);
        a0 += shift0;
        rsum -= shift0 * n;
        for (int k = 0; k < K; ++k)
          for (int d = 0; d < D; ++d) H[d + D * k] -= shift0 * S[d + D * k];

        for (int g = 0; g < J; ++g) {
          if (!active[g]) continue;
          fit_group(g, ell, v, max_change);
        }

        a = bcoef;
        if (max_change < eps) {
          active_converged = true;
          break;
        }
      }

      if (!active_converged) break;

      // Scan the other SNPs and add the nonzero ones to the active set
      int violations = 0;
      double max_change_scan = 0.0;
      for (int g = 0; g < J; ++g) {
        if (active[g]) continue;
        if (fit_group(g, ell, v, max_change_scan)) {
          active[g] = 1;
          violations++;
        }
      }

      a = bcoef;
      if (violations == 0) {
        converged[ell] = true;
        break;
      }
    }

    intercept[ell] = a0;
    for (int j = 0; j < p; ++j) beta(j, ell) = bcoef[j];
  }

  return List::create(
      _["beta_ortho"] = beta,
      _["intercept_ortho"] = intercept,
      _["iter"] = iter,
      _["converged"] = converged);
}
