#define TMB_LIB_INIT R_init_jointCPUE
#define EIGEN_DONT_PARALLELIZE
#include <TMB.hpp>
#define _USE_MATH_DEFINES
#include <cmath>

template<class Type>
Type dlnorm_bc(const Type& x, const Type& meanlog, const Type& sdlog, int give_log = 0) {
  Type adjusted_meanlog = meanlog - pow(sdlog, 2) / 2;
  Type logres = dnorm(log(x), adjusted_meanlog, sdlog, true) - log(x);

  if (give_log) {
    return logres;
  } else {
    return exp(logres);
  }
}

template <class Type>
Type pc_prior_matern(
  Type range,
  Type sigma,
  Type matern_range,
  Type matern_sigma,
  Type range_prob,
  Type sigma_prob,
  int give_log = 0,
  int share_range = 0
) {
  Type d = 2.;
  Type dhalf = d / 2.;

  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(sigma_prob) / matern_sigma;

  Type range_ll = log(dhalf)
    + log(lam1)
    + (-1.0 - dhalf) * log(range)
    - lam1 * pow(range, -dhalf);

  Type sigma_ll = log(lam2) - lam2 * sigma;

  Type penalty = sigma_ll;
  if (!share_range) penalty += range_ll;

  if (give_log) {
    return penalty;
  } else {
    return exp(penalty);
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_INTEGER(n_i);
  DATA_INTEGER(n_t);
  DATA_INTEGER(n_f);
  DATA_INTEGER(n_g);
  DATA_INTEGER(n_m);

  DATA_VECTOR(b_i);
  DATA_IVECTOR(t_i);
  DATA_IVECTOR(f_i);
  DATA_IVECTOR(month_i);
  DATA_IMATRIX(has_tf);
  DATA_VECTOR(area_g);

  DATA_SPARSE_MATRIX(A_is);
  DATA_SPARSE_MATRIX(A_gs);
  DATA_IMATRIX(Ais_ij);
  DATA_VECTOR(Ais_x);

  DATA_SCALAR(matern_range);
  DATA_SCALAR(range_prob);
  DATA_SCALAR(matern_sigma_0);
  DATA_SCALAR(matern_sigma_t);
  DATA_SCALAR(matern_sigma_fleet);
  DATA_SCALAR(sigma_prob);
  DATA_INTEGER(use_pop_spatial);
  DATA_INTEGER(use_pop_spatiotemporal);
  DATA_INTEGER(use_pop_spatiotemporal_rw);
  DATA_INTEGER(use_month_fe);
  DATA_INTEGER(use_fleet_sd);
  DATA_INTEGER(use_q_diffs_time);
  DATA_INTEGER(use_q_diffs_spatial);

  DATA_STRUCT(spde, spde_aniso_t);

  PARAMETER(ln_sd);
  PARAMETER_VECTOR(ln_sd_fleet);
  PARAMETER_VECTOR(ln_H_input);

  PARAMETER(ln_range_1);
  PARAMETER(ln_sigma_0_1);
  PARAMETER(ln_sigma_t_1);

  PARAMETER_VECTOR(yq_t_1);
  PARAMETER_VECTOR(month_beta);
  PARAMETER_VECTOR(omega_s_1);
  PARAMETER_MATRIX(epsilon_st_1);
  PARAMETER_VECTOR(fleet_f);
  PARAMETER(fleet_ln_std_dev);
  PARAMETER_MATRIX(fleet_t);
  PARAMETER(fleet_t_ln_std_dev);
  PARAMETER_MATRIX(fleet_s);
  PARAMETER(ln_sigma_fleet);
  PARAMETER_VECTOR(eps_index);

  Type nll = 0;
  Type nll_prior = 0;
  Type nll_penalty = 0;
  Type sd_shared = exp(ln_sd);
  Type range_1 = exp(ln_range_1);
  Type sigma_0_1 = exp(ln_sigma_0_1);
  Type sigma_t_1 = exp(ln_sigma_t_1);
  Type fleet_std_dev = exp(fleet_ln_std_dev);
  Type fleet_t_std_dev = exp(fleet_t_ln_std_dev);
  vector<Type> sd_fleet(n_f);
  bool range_prior_added = false;
  Type kappa_1 = sqrt(Type(8.0)) / range_1;

  for (int f = 0; f < n_f; f++) {
    sd_fleet(f) = (use_fleet_sd == 1) ? exp(ln_sd_fleet(f)) : sd_shared;
  }

  matrix<Type> H(2, 2);
  H(0, 0) = exp(ln_H_input(0));
  H(1, 0) = ln_H_input(1);
  H(0, 1) = ln_H_input(1);
  H(1, 1) = (1 + ln_H_input(1) * ln_H_input(1)) / exp(ln_H_input(0));

  SparseMatrix<Type> Q_1 = Q_spde(spde, kappa_1, H);

  if (use_pop_spatial == 1) {
    nll_prior -= pc_prior_matern(
      range_1, sigma_0_1,
      matern_range, matern_sigma_0,
      range_prob, sigma_prob,
      1, range_prior_added ? 1 : 0
    );
    range_prior_added = true;

    Type tau_0_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_0_1);
    nll += SCALE(GMRF(Q_1), 1. / tau_0_1)(omega_s_1);
  }

  if (use_pop_spatiotemporal == 1) {
    nll_prior -= pc_prior_matern(
      range_1, sigma_t_1,
      matern_range, matern_sigma_t,
      range_prob, sigma_prob,
      1, range_prior_added ? 1 : 0
    );
    range_prior_added = true;

    Type tau_t_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_t_1);
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_1), 1. / tau_t_1)(epsilon_st_1.col(0));
      for (int t = 1; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), 1. / tau_t_1)(epsilon_st_1.col(t) - epsilon_st_1.col(t - 1));
      }
    } else {
      for (int t = 0; t < n_t; t++) {
        nll += SCALE(GMRF(Q_1), 1. / tau_t_1)(epsilon_st_1.col(t));
      }
    }
  }

  for (int f = 0; f < n_f - 1; f++) {
    nll -= dnorm(fleet_f(f), Type(0.0), fleet_std_dev, true);
  }

  if (use_q_diffs_time == 1 && n_f > 1) {
    for (int j = 0; j < n_f - 1; j++) {
      for (int t = 0; t < n_t; t++) {
        if (has_tf(t, j) == 1) {
          nll -= dnorm(fleet_t(t, j), Type(0.0), fleet_t_std_dev, true);
        }
      }
    }
  }

  if (use_q_diffs_spatial == 1 && n_f > 1) {
    Type sigma_fleet = exp(ln_sigma_fleet);
    Type tau_fleet = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_fleet);

    nll_prior -= pc_prior_matern(
      range_1, sigma_fleet,
      matern_range, matern_sigma_fleet,
      range_prob, sigma_prob,
      1, range_prior_added ? 1 : 0
    );
    range_prior_added = true;

    for (int f = 0; f < n_f - 1; f++) {
      nll += SCALE(GMRF(Q_1), 1. / tau_fleet)(fleet_s.col(f));
    }
  }

  vector<Type> s_effect_1(n_i);
  vector<Type> st_effect_1(n_i);
  vector<Type> month_effect_m(n_m);
  vector<Type> fleet_t_mean(n_f > 1 ? n_f - 1 : 0);
  vector<Type> fleet_s_effect(n_i);
  vector<Type> eta_hat_i(n_i);
  vector<Type> mu_hat_i(n_i);
  vector<Type> residual_raw_i(n_i);
  vector<Type> residual_log_i(n_i);
  s_effect_1.setZero();
  st_effect_1.setZero();
  month_effect_m.setZero();
  fleet_t_mean.setZero();
  fleet_s_effect.setZero();
  eta_hat_i.setZero();
  mu_hat_i.setZero();
  residual_raw_i.setZero();
  residual_log_i.setZero();

  if (use_pop_spatial == 1) {
    s_effect_1 = A_is * omega_s_1;
  }
  
  if (use_month_fe == 1 && n_m > 1) {
    Type sum_month = 0.0;
    for (int m = 0; m < n_m - 1; m++) {
      month_effect_m(m) = month_beta(m);
      sum_month += month_beta(m);
    }
    month_effect_m(n_m - 1) = -sum_month;
  }

  if (use_q_diffs_time == 1 && n_f > 1) {
    for (int j = 0; j < n_f - 1; j++) {
      Type sum = 0.0;
      Type cnt = 0.0;
      for (int t = 0; t < n_t; t++) {
        if (has_tf(t, j) == 1) {
          sum += fleet_t(t, j);
          cnt += 1.0;
        }
      }
      if (cnt > 0.0) {
        fleet_t_mean(j) = sum / cnt;
      }
    }
  }

  for (int r = 0; r < Ais_ij.rows(); r++) {
    int i = Ais_ij(r, 0);
    int s = Ais_ij(r, 1);
    int t_id = t_i(i);
    int f_id = f_i(i);
    if (use_pop_spatiotemporal == 1) {
      st_effect_1(i) += Ais_x(r) * epsilon_st_1(s, t_id);
    }
    if (use_q_diffs_spatial == 1 && f_id > 0) {
      fleet_s_effect(i) += Ais_x(r) * fleet_s(s, f_id - 1);
    }
  }

  for (int i = 0; i < n_i; i++) {
    int tid = t_i(i);
    int fid = f_i(i);

    Type yq_effect_1 = yq_t_1(tid);
    Type month_effect = Type(0.0);
    if (use_month_fe == 1 && n_m > 1) {
      month_effect = month_effect_m(month_i(i));
    }
    Type f_effect_1 = (fid == 0) ? Type(0) : fleet_f(fid - 1);
    Type fleet_t_effect = Type(0.0);
    if (use_q_diffs_time == 1 && fid > 0) {
      int j = fid - 1;
      if (has_tf(tid, j) == 1) {
        fleet_t_effect = fleet_t(tid, j) - fleet_t_mean(j);
      }
    }
    Type eta1 = yq_effect_1 + month_effect + s_effect_1(i) + st_effect_1(i) + f_effect_1 + fleet_t_effect + fleet_s_effect(i);
    eta_hat_i(i) = eta1;
    mu_hat_i(i) = exp(eta1);
    residual_raw_i(i) = b_i(i) - mu_hat_i(i);
    residual_log_i(i) = log(b_i(i)) - eta1;
    nll -= dlnorm_bc(b_i(i), eta1, sd_fleet(fid), true);
  }

  vector<Type> s_effect_proj_1(n_g);
  matrix<Type> st_effect_proj_1(n_g, n_t);
  s_effect_proj_1.setZero();
  st_effect_proj_1.setZero();

  if (use_pop_spatial == 1) {
    s_effect_proj_1 = A_gs * omega_s_1;
  }
  if (use_pop_spatiotemporal == 1) {
    st_effect_proj_1 = A_gs * epsilon_st_1;
  }

  // Standardized index: sum over grid of (area * exp(eta)) per time t.
  // When use_month_fe == 1 (yearly index), month fixed effect is deliberately
  // excluded here so the index reflects year-level abundance only.
  matrix<Type> cpue_density(n_g, n_t);
  vector<Type> mu_total(n_t);
  vector<Type> link_total(n_t);
  mu_total.setZero();
  Type max_eta_proj = Type(40.0);

  for (int t = 0; t < n_t; t++) {
    Type yq_effect_proj_1 = yq_t_1(t);
    for (int g = 0; g < n_g; g++) {
      Type eta1_proj = yq_effect_proj_1 + s_effect_proj_1(g) + st_effect_proj_1(g, t);
      Type safe_eta1_proj = CppAD::CondExpGt(eta1_proj, max_eta_proj, max_eta_proj, eta1_proj);
      safe_eta1_proj = CppAD::CondExpLt(safe_eta1_proj, -max_eta_proj, -max_eta_proj, safe_eta1_proj);
      Type cpue = exp(safe_eta1_proj);
      cpue_density(g, t) = cpue;
      mu_total(t) += cpue * area_g(g);
    }
    link_total(t) = log(mu_total(t) + Type(1e-12));
  }

  if (eps_index.size() > 0) {
    Type S;
    for (int t = 0; t < n_t; t++) {
      S = mu_total(t);
      S = newton::Tag(S);
      nll_penalty += eps_index(t) * S;
    }
  }

  nll += nll_prior;
  nll += nll_penalty;

  REPORT(cpue_density);
  REPORT(fleet_f);
  REPORT(fleet_t);
  REPORT(fleet_s);
  REPORT(month_effect_m);
  REPORT(sd_fleet);
  REPORT(eta_hat_i);
  REPORT(mu_hat_i);
  REPORT(residual_raw_i);
  REPORT(residual_log_i);
  ADREPORT(link_total);
  REPORT(nll_prior);
  REPORT(nll_penalty);
  return nll;
}
