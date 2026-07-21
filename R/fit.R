#' @useDynLib jointCPUE, .registration = TRUE
NULL

#' Fit an integrated spatiotemporal CPUE standardization model
#'
#' Fits an integrated spatiotemporal model for CPUE standardization that jointly
#' models observation and sampling processes across one or more fisheries or surveys.
#' The framework is designed to reduce bias caused by preferential sampling,
#' targeting behavior, and heterogeneous effort distributions, and to estimate
#' a coherent latent abundance index.
#'
#' The model is implemented using Template Model Builder (TMB). Spatial and
#' spatiotemporal random fields are represented using the SPDE (stochastic partial
#' differential equation) approach for computational efficiency. See the package
#' vignettes for worked examples and model-comparison workflows.
#'
#' @param data_utm A data.frame with required columns.
#' @param mesh An `jointCPUEmesh` or a bare fmesher mesh.
#' @param formula_catchability Optional one-sided formula defining observation
#'   catchability covariates. `s()` terms are penalized smooths; simple numeric
#'   terms are linear fixed effects; simple factor/character terms are
#'   sum-to-zero fixed effects; `f(var)` terms are random factor effects. These
#'   terms enter the observation model only. Numeric category codes must be
#'   converted to factors to request a fixed factor effect.
#' @param formula_population Optional one-sided formula defining population
#'   covariates. These terms enter both the observation model and the projected
#'   standardized index.
#' @param projection_data Optional extrapolation-grid covariates for
#'   `formula_population`. Supply one row per grid cell for static covariates,
#'   or one row per grid cell-time combination with `tid` for time-varying
#'   covariates.
#' @param pop_spatial "on" or "off". If "off", the population-level
#'   time-constant spatial field is excluded from both the likelihood and index.
#' @param pop_spatiotemporal "on" or "off". If "off", the population-level
#'   spatiotemporal field is excluded from both the likelihood and index.
#' @param pop_spatiotemporal_type Either `"iid"` or `"rw"`. Controls the
#'   temporal structure of the population-level spatiotemporal field when
#'   `pop_spatiotemporal = "on"`. If `"iid"`, each time slice is independent.
#'   If `"rw"`, the field follows a first-order random walk through time.
#' @param q_diffs_system "on" or "off". If "off", fleet effects are fixed to zero.
#' @param q_diffs_time "on" or "off". If "on", estimate fleet-specific temporal deviations.
#' @param q_diffs_spatial "on" or "off". If "on", estimate fleet-specific spatial deviations.
#' @param index Either `"monthly"` or `"yearly"`. Chooses the time resolution of the
#'   standardized index.
#' @param month_col Optional column name for calendar month (`1`--`12`). Used only
#'   when `index = "yearly"` and `month_diffs = "on"`. If `NULL`, `jointCPUE()`
#'   looks for `month` or `Month`.
#' @param month_diffs "on" or "off". Controls whether a calendar month fixed
#'   effect is added when `index = "yearly"`. The month effect uses a sum-to-zero
#'   constraint over observed month levels, enters the observation model only,
#'   and is excluded from the standardized index. When `index = "monthly"`, the
#'   month effect is always disabled regardless of this setting.
#' @param obs_sd Either `"shared"` or `"fleet"`. Controls the observation-error
#'   standard deviation in the lognormal likelihood. `"shared"` uses a single
#'   SD across all fleets. `"fleet"` estimates one SD per fleet.
#' @param extrapolation_grid Optional data.frame defining a common
#'   extrapolation grid/key across fits. Use this when comparing joint and
#'   single-fleet models on the same spatial integration domain. The supplied
#'   object must contain `utm_x_scale` and `utm_y_scale`, plus either `area_km2`
#'   or enough information to compute area (`lon`/`lat`, or `lon_std`/`lat`).
#' @param control Control list passed to [stats::nlminb()].
#' @param ncores Optional integer. If provided, sets the number of OpenMP threads. Passed to [TMB::openmp()].
#' @param silent Logical. Passed to [TMB::MakeADFun()].
#' @param ... Passed to [jointCPUE::make_data()] (for example `area_scale`).
#'
#' @return An object of class `jointCPUE` with elements `obj`, `opt`, `rep`, `prep`, etc.
#' @author Rujia Bi \email{bikayla5@gmail.com}
#' @export
jointCPUE <- function(
    data_utm,
    mesh,
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    pop_spatial = c("on", "off"),
    pop_spatiotemporal = c("on", "off"),
    pop_spatiotemporal_type = c("iid", "rw"),
    q_diffs_system = c("on", "off"),
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    control = list(eval.max = 1e5, iter.max = 1e5),
    ncores = NULL,
    ...,
    index = c("monthly", "yearly"),
    month_col = NULL,
    month_diffs = c("on", "off"),
    obs_sd = c("shared", "fleet"),
    extrapolation_grid = NULL,
    silent = FALSE
) {
  pop_spatial <- match.arg(pop_spatial)
  pop_spatiotemporal <- match.arg(pop_spatiotemporal)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  q_diffs_system  <- match.arg(q_diffs_system)
  q_diffs_time    <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  index <- match.arg(index)
  month_diffs <- match.arg(month_diffs)
  obs_sd <- match.arg(obs_sd)
  month_diffs <- if (index == "yearly") month_diffs else "off"
  
  data_utm <- as.data.frame(data_utm)
  
  # ---- 1) Data prep (single source of truth) ----
  prep <- make_data(
    data_utm = data_utm,
    mesh = mesh,
    index = index,
    month_col = month_col,
    month_diffs = month_diffs,
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    projection_data = projection_data,
    extrapolation_grid = extrapolation_grid,
    ...
  )
  
  data_tmb <- prep$data
  
  # ---- 2) Defensive checks (catch mismatches early) ----
  n_t <- data_tmb$n_t
  n_f <- data_tmb$n_f
  n_s <- data_tmb$spde$n_s
  data_tmb$use_pop_spatial <- as.integer(pop_spatial == "on")
  data_tmb$use_pop_spatiotemporal <- as.integer(pop_spatiotemporal == "on")
  data_tmb$use_pop_spatiotemporal_rw <- as.integer(
    pop_spatiotemporal == "on" && pop_spatiotemporal_type == "rw"
  )
  data_tmb$use_month_fe <- as.integer(
    index == "yearly" && month_diffs == "on" && data_tmb$n_m > 1L
  )
  data_tmb$use_fleet_sd <- as.integer(obs_sd == "fleet")
  data_tmb$use_q_diffs_system <- as.integer(q_diffs_system == "on" && n_f > 1L)
  data_tmb$use_q_diffs_time <- as.integer(q_diffs_time == "on" && n_f > 1L)
  data_tmb$use_q_diffs_spatial <- as.integer(q_diffs_spatial == "on" && n_f > 1L)

  has_smooths_catch <- isTRUE(data_tmb$has_smooths_catch == 1L)
  K_smooth_catch <- if (has_smooths_catch) ncol(data_tmb$Xs_catch) else 0L
  n_smooth_catch <- if (has_smooths_catch) length(data_tmb$Zs_catch) else 0L
  sum_k_catch <- if (has_smooths_catch && n_smooth_catch > 0L) {
    sum(vapply(data_tmb$Zs_catch, ncol, 0L))
  } else {
    0L
  }
  K_fixed_catch <- if (isTRUE(data_tmb$has_fixed_catch == 1L)) ncol(data_tmb$Xf_catch) else 0L
  has_random_factors_catch <- isTRUE(data_tmb$has_random_factors_catch == 1L)
  n_factor_random_catch <- if (has_random_factors_catch) length(data_tmb$Zf_catch) else 0L
  sum_k_factor_catch <- if (has_random_factors_catch && n_factor_random_catch > 0L) {
    sum(vapply(data_tmb$Zf_catch, ncol, 0L))
  } else {
    0L
  }

  has_smooths_pop <- isTRUE(data_tmb$has_smooths_pop == 1L)
  K_smooth_pop <- if (has_smooths_pop) ncol(data_tmb$Xs_pop_i) else 0L
  n_smooth_pop <- if (has_smooths_pop) length(data_tmb$Zs_pop_i) else 0L
  sum_k_pop <- if (has_smooths_pop && n_smooth_pop > 0L) {
    sum(vapply(data_tmb$Zs_pop_i, ncol, 0L))
  } else {
    0L
  }
  K_fixed_pop <- if (isTRUE(data_tmb$has_fixed_pop == 1L)) ncol(data_tmb$Xf_pop_i) else 0L
  has_random_factors_pop <- isTRUE(data_tmb$has_random_factors_pop == 1L)
  n_factor_random_pop <- if (has_random_factors_pop) length(data_tmb$Zf_pop_i) else 0L
  sum_k_factor_pop <- if (has_random_factors_pop && n_factor_random_pop > 0L) {
    sum(vapply(data_tmb$Zf_pop_i, ncol, 0L))
  } else {
    0L
  }
  
  # ---- 3) Initial parameters (must match cpp) ----
  parameters <- .make_parameters_jointCPUE(
    n_t = n_t,
    n_f = n_f,
    n_s = n_s,
    n_m = data_tmb$n_m,
    K_smooth_catch = K_smooth_catch,
    n_smooth_catch = n_smooth_catch,
    sum_k_catch = sum_k_catch,
    K_smooth_pop = K_smooth_pop,
    n_smooth_pop = n_smooth_pop,
    sum_k_pop = sum_k_pop,
    K_fixed_catch = K_fixed_catch,
    K_fixed_pop = K_fixed_pop,
    n_factor_random_catch = n_factor_random_catch,
    sum_k_factor_catch = sum_k_factor_catch,
    n_factor_random_pop = n_factor_random_pop,
    sum_k_factor_pop = sum_k_factor_pop
  )
  
  # ---- 4) MAP (turn on/off components without touching cpp) ----
  map <- .make_map_jointCPUE(
    parameters = parameters,
    n_f = n_f,
    pop_spatial = pop_spatial,
    pop_spatiotemporal = pop_spatiotemporal,
    q_diffs_system = q_diffs_system,
    q_diffs_time = q_diffs_time,
    q_diffs_spatial = q_diffs_spatial,
    index = index,
    month_diffs = month_diffs,
    obs_sd = obs_sd,
    n_m = data_tmb$n_m,
    has_tf = data_tmb$has_tf
  )
  
  # ---- 5) Random effects list ----
  random <- character(0)

  if (pop_spatial == "on") {
    random <- c(random, "omega_s_1")
  }

  if (pop_spatiotemporal == "on") {
    random <- c(random, "epsilon_st_1")
  }

  if (n_f > 1L && q_diffs_system == "on") {
    random <- c(random, "fleet_f")
  }

  if (n_f > 1L && q_diffs_time == "on" && any(data_tmb$has_tf > 0L)) {
    random <- c(random, "fleet_t")
  }

  if (n_f > 1L && q_diffs_spatial == "on") {
    random <- c(random, "fleet_s")
  }

  if (sum_k_catch > 0L) {
    random <- c(random, "b_smooth_catch")
  }
  if (sum_k_pop > 0L) {
    random <- c(random, "b_smooth_pop")
  }
  if (sum_k_factor_catch > 0L) {
    random <- c(random, "b_factor_catch")
  }
  if (sum_k_factor_pop > 0L) {
    random <- c(random, "b_factor_pop")
  }
  
  random <- unique(random)
  
  # ---- 6) Build & optimize ----
  .check_fit_inputs_jointCPUE(data_tmb)
  
  DLL <- "jointCPUE"
  
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      stop("`ncores` must be a positive integer.", call. = FALSE)
    }
    if (ncores > 1L) {
      TMB::openmp(ncores, autopar = TRUE, DLL = DLL)
    }
  }
  
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    map = map,
    random = random,
    DLL = DLL,
    silent = silent
  )
  
  opt <- .safe_optimize(obj, control)
  
  rep <- TMB::sdreport(obj)
  
  out <- list(
    obj = obj,
    opt = opt,
    rep = rep,
    prep = prep,
    data_tmb = data_tmb,
    map = map,
    random = random,
    control = control,
    settings = list(
      pop_spatial = pop_spatial,
      pop_spatiotemporal = pop_spatiotemporal,
      pop_spatiotemporal_type = pop_spatiotemporal_type,
      q_diffs_system = q_diffs_system,
      q_diffs_time = q_diffs_time,
      q_diffs_spatial = q_diffs_spatial,
      formula_catchability = formula_catchability,
      formula_population = formula_population,
      index = index,
      month_diffs = month_diffs,
      obs_sd = obs_sd,
      extrapolation_grid_supplied = !is.null(extrapolation_grid),
      month_col = prep$time$month_col,
      month_values = prep$time$month_values,
      DLL = DLL,
      ncores = ncores
    ),
    diagnostics = list(
      convergence = opt$convergence,
      message = opt$message,
      max_grad = max(abs(obj$gr(opt$par)))
    )
  )
  class(out) <- "jointCPUE"
  out
}
