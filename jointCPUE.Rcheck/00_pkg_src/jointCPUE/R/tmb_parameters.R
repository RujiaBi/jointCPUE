#' Internal: Build initial parameter list for the jointCPUE TMB model
#'
#' @param n_t Number of time steps.
#' @param n_f Number of fleets including the baseline fleet coded as 0.
#' @param n_s Number of SPDE mesh vertices.
#' @param n_m Number of observed calendar month levels used by the yearly
#'   month fixed effect.
#'
#' @return Named list of initial parameter values for TMB::MakeADFun().
#' @keywords internal
.make_parameters_jointCPUE <- function(n_t, n_f, n_s, n_m = 1L) {
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
