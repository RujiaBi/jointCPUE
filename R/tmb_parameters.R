#' Internal: Build initial parameter list for the jointCPUE TMB model
#'
#' @param n_t Number of time steps.
#' @param n_f Number of fleets including the baseline fleet coded as 0.
#' @param n_s Number of SPDE mesh vertices.
#' @param n_m Number of observed calendar month levels used by the yearly
#'   month fixed effect.
#' @param K_smooth_catch Number of fixed-effect columns from catchability
#'   smooth terms.
#' @param n_smooth_catch Number of catchability smooth penalty blocks.
#' @param sum_k_catch Total penalized coefficients across catchability smooths.
#' @param K_smooth_pop Number of fixed-effect columns from population smooths.
#' @param n_smooth_pop Number of population smooth penalty blocks.
#' @param sum_k_pop Total penalized coefficients across population smooths.
#' @param K_fixed_catch Number of catchability fixed covariate columns.
#' @param K_fixed_pop Number of population fixed covariate columns.
#' @param n_factor_random_catch Number of catchability random-factor terms.
#' @param sum_k_factor_catch Total random-factor coefficients for catchability.
#' @param n_factor_random_pop Number of population random-factor terms.
#' @param sum_k_factor_pop Total random-factor coefficients for population.
#'
#' @return Named list of initial parameter values for TMB::MakeADFun().
#' @keywords internal
.make_parameters_jointCPUE <- function(
    n_t,
    n_f,
    n_s,
    n_m = 1L,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L,
    K_fixed_catch = 0L,
    K_fixed_pop = 0L,
    n_factor_random_catch = 0L,
    sum_k_factor_catch = 0L,
    n_factor_random_pop = 0L,
    sum_k_factor_pop = 0L
) {
  empty_mat <- function(nr, nc) {
    matrix(0.0, nrow = nr, ncol = nc)
  }

  list(
    ln_sd = 0.0,
    ln_sd_fleet = rep(0.0, n_f),
    ln_H_input = c(0.0, 0.0),
    ln_range_1 = 0.0,
    ln_sigma_0_1 = 0.0,
    ln_sigma_t_1 = 0.0,
    yq_t_1 = rep(0.0, n_t),
    month_beta = if (n_m > 1L) rep(0.0, n_m - 1L) else numeric(0),
    bs_catch = rep(0.0, K_smooth_catch),
    b_smooth_catch = rep(0.0, sum_k_catch),
    ln_smooth_sigma_catch = rep(0.0, n_smooth_catch),
    bf_catch = rep(0.0, K_fixed_catch),
    b_factor_catch = rep(0.0, sum_k_factor_catch),
    ln_factor_sigma_catch = rep(0.0, n_factor_random_catch),
    bs_pop = rep(0.0, K_smooth_pop),
    b_smooth_pop = rep(0.0, sum_k_pop),
    ln_smooth_sigma_pop = rep(0.0, n_smooth_pop),
    bf_pop = rep(0.0, K_fixed_pop),
    b_factor_pop = rep(0.0, sum_k_factor_pop),
    ln_factor_sigma_pop = rep(0.0, n_factor_random_pop),
    omega_s_1 = rep(0.0, n_s),
    epsilon_st_1 = empty_mat(n_s, n_t),
    fleet_f = if (n_f > 1L) rep(0.0, n_f - 1L) else numeric(0),
    fleet_ln_std_dev = 0.0,
    fleet_t = if (n_f > 1L) empty_mat(n_t, n_f - 1L) else empty_mat(n_t, 0L),
    fleet_t_ln_std_dev = 0.0,
    fleet_s = if (n_f > 1L) empty_mat(n_s, n_f - 1L) else empty_mat(n_s, 0L),
    ln_sigma_fleet = 0.0,
    eps_index = numeric(0)
  )
}
