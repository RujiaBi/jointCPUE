#' Prepare data objects and mesh for jointCPUE workflows
#'
#' Main data-prep function for jointCPUE.
#' - mesh/SPDE/A matrices
#' - extrapolation key grid + areas
#'
#' @param data_utm A data.frame containing required columns (with utm_x/y_scale)
#' @param mesh jointCPUEmesh built from make_mesh(), or a custom mesh
#' @param area_scale Numeric or "auto". Scaling factor for area_km2.
#' @param index Either `"monthly"` or `"yearly"`.
#' @param month_col Optional column name for calendar month (`1`--`12`). Used
#'   only when `index = "yearly"` and `month_diffs = "on"`. If `NULL`, looks for
#'   `month` or `Month`.
#' @param month_diffs "on" or "off". Controls whether a calendar month fixed
#'   effect is prepared when `index = "yearly"`. The month effect uses a
#'   sum-to-zero constraint over observed month levels, enters the observation
#'   model only, and is excluded from the standardized index. When
#'   `index = "monthly"`, the month effect is always disabled.
#' @param formula_catchability Optional one-sided formula defining observation
#'   catchability covariates. `s()` terms are treated as penalized smooths;
#'   simple numeric terms are linear fixed effects; simple factor/character
#'   terms are fixed effects with sum-to-zero coding; `f(var)` terms are random
#'   factor effects. These terms enter the observation model only and are
#'   excluded from the standardized index. Numeric category codes must be
#'   converted to factors to request a fixed factor effect.
#' @param formula_population Optional one-sided formula defining population
#'   covariates. These terms enter both the observation model and the projected
#'   standardized index.
#' @param projection_data Optional data.frame giving extrapolation-grid
#'   covariates for `formula_population`. Static covariates require one row per
#'   grid cell with `utm_x_scale`, `utm_y_scale`, and covariate columns.
#'   Time-varying covariates require one row per grid cell-time combination and
#'   a `tid` column using the same 0-based coding as `data_utm$tid`.
#' @param extrapolation_grid Optional data.frame defining the extrapolation
#'   grid/key shared across fits. Use this to force joint and single-fleet
#'   models onto the same spatial integration domain. If `NULL`, the grid is
#'   built automatically from unique observed locations in `data_utm`. Supplied
#'   grids must contain `utm_x_scale` and `utm_y_scale`, plus either `area_km2`
#'   or enough information to compute area (`lon`/`lat`, or `lon_std`/`lat`).
#'
#' @return A list with elements `data`, `key`, `scales`, and `time`.
#' @author Rujia Bi \email{bikayla5@gmail.com}
#' @export
make_data <- function(
    data_utm,
    mesh,
    area_scale = "auto",
    index = c("monthly", "yearly"),
    month_col = NULL,
    month_diffs = c("on", "off"),
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    extrapolation_grid = NULL
) {
  index <- match.arg(index)
  month_diffs <- match.arg(month_diffs)
  month_diffs <- if (index == "yearly") month_diffs else "off"
  data_utm <- as.data.frame(data_utm)
  
  .check_required_cols(data_utm, c("cpue", "lon", "lat", "tid", "fleetid", "utm_x_scale", "utm_y_scale"))
  .check_numeric(data_utm, c("cpue", "lon", "lat", "tid", "fleetid", "utm_x_scale", "utm_y_scale"))
  
  if (anyNA(data_utm$lon) || anyNA(data_utm$lat)) {
    stop("`lon`/`lat` must not contain NA.", call. = FALSE)
  }
  
  if (anyNA(data_utm$utm_x_scale) || anyNA(data_utm$utm_y_scale)) {
    stop("`utm_x_scale`/`utm_y_scale` must not contain NA.", call. = FALSE)
  }

  if (any(!is.finite(data_utm$cpue)) || any(data_utm$cpue <= 0)) {
    stop(
      "`cpue` must be strictly positive and finite. Add a small floor to zero values before fitting.",
      call. = FALSE
    )
  }
  
  # ---- SPDE + A matrix (handle jointCPUEmesh or bare mesh) ----
  loc_xy <- as.matrix(data_utm[, c("utm_x_scale", "utm_y_scale"), drop = FALSE])
  
  # mesh can be:
  #  - jointCPUEmesh from make_mesh()
  #  - a bare fmesher mesh object (mesh$mesh)
  mesh_in <- mesh
  mesh_obj <- .as_jointCPUEmesh(
    mesh = mesh_in,
    loc_xy = loc_xy,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    recompute_A = "auto"
  )
  
  if (!is.null(mesh_obj$loc_xy)) {
    r1 <- range(mesh_obj$loc_xy[,1]); r2 <- range(loc_xy[,1])
    if (is.finite(r1[1]) && is.finite(r2[1])) {
      if (abs(diff(r1) - diff(r2)) / max(1e-12, diff(r2)) > 0.5) {
        warning("Mesh coordinate scale may not match `utm_x_scale/utm_y_scale`. Check scaling.", call. = FALSE)
      }
    }
  }
  
  mesh <- mesh_obj$mesh
  
  # ---- SPDE (INLA) ----
  spde <- mesh_obj$spde
  
  # ---- A matrices ----
  A_is <- mesh_obj$A
  A_isT <- methods::as(A_is, "TsparseMatrix")
  Ais_ij <- cbind(A_isT@i, A_isT@j)
  Ais_x  <- A_isT@x
  
  # ---- key/extrapolation grid ----
  key_out <- .prep_key_area(
    data_utm,
    mesh,
    area_scale = area_scale,
    extrapolation_grid = extrapolation_grid
  )
  key <- key_out$key
  A_gs <- key_out$A_gs
  
  n_i <- nrow(data_utm)
  
  # ---- user-supplied 0-based indices: validate only ----
  t_chk <- .check_0based_contiguous(data_utm$tid, "tid")
  f_chk <- .check_0based_contiguous(data_utm$fleetid, "fleetid")
  month_prep <- .prep_month_index(
    data_utm,
    index = index,
    month_col = month_col,
    month_diffs = month_diffs
  )
  .check_month_formula_conflict(
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    month_col = month_prep$month_col,
    use_month_fe = month_prep$use_month_fe
  )
  
  t_i <- t_chk$x
  f_i <- f_chk$x
  
  n_t <- t_chk$n
  n_f <- f_chk$n
  tid_values <- seq.int(0L, n_t - 1L)

  has_tf <- matrix(0L, nrow = n_t, ncol = max(0L, n_f - 1L))
  if (n_f > 1L) {
    ii <- which(f_i > 0L)
    if (length(ii) > 0L) {
      tt <- t_i[ii] + 1L
      ff <- f_i[ii]
      has_tf[cbind(tt, ff)] <- 1L
    }
  }

  sm_catch <- .normalize_smoother_output(
    parse_smoothers(
      formula = formula_catchability,
      data = data_utm,
      knots = NULL,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )
  fx_catch <- .normalize_fixed_output(
    parse_fixed_covariates(
      formula = formula_catchability,
      data = data_utm,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )
  rf_catch <- .normalize_random_factor_output(
    parse_random_factors(
      formula = formula_catchability,
      data = data_utm,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )

  sm_pop_obs <- .normalize_smoother_output(
    parse_smoothers(
      formula = formula_population,
      data = data_utm,
      knots = NULL,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )
  fx_pop_obs <- .normalize_fixed_output(
    parse_fixed_covariates(
      formula = formula_population,
      data = data_utm,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )
  rf_pop_obs <- .normalize_random_factor_output(
    parse_random_factors(
      formula = formula_population,
      data = data_utm,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )

  if (isTRUE(sm_pop_obs$has_smooths) || isTRUE(fx_pop_obs$has_fixed) || isTRUE(rf_pop_obs$has_random)) {
    needed_vars_pop <- unique(c(
      .smooth_vars_from_basis(sm_pop_obs$basis_out),
      .fixed_vars_from_basis(fx_pop_obs$basis_out),
      .random_factor_vars_from_basis(rf_pop_obs$basis_out)
    ))
    projection_data_use <- if (is.null(projection_data)) {
      .default_population_projection_data(
        data_utm = data_utm,
        key = key,
        needed_vars = needed_vars_pop,
        tid_values = tid_values
      )
    } else {
      .align_projection_data_to_key(
        projection_data = projection_data,
        key = key,
        needed_vars = needed_vars_pop,
        tid_values = tid_values
      )
    }
    .warn_projection_na(projection_data_use, needed_vars_pop)

    sm_pop_proj <- .normalize_smoother_output(
      parse_smoothers(
        formula = formula_population,
        data = data_utm,
        knots = NULL,
        newdata = projection_data_use,
        basis_prev = sm_pop_obs$basis_out
      ),
      n_rows = nrow(key) * n_t
    )
    fx_pop_proj <- .normalize_fixed_output(
      parse_fixed_covariates(
        formula = formula_population,
        data = data_utm,
        newdata = projection_data_use,
        basis_prev = fx_pop_obs$basis_out
      ),
      n_rows = nrow(key) * n_t
    )
    rf_pop_proj <- .normalize_random_factor_output(
      parse_random_factors(
        formula = formula_population,
        data = data_utm,
        newdata = projection_data_use,
        basis_prev = rf_pop_obs$basis_out
      ),
      n_rows = nrow(key) * n_t
    )
  } else {
    projection_data_use <- .expand_projection_over_time(
      key[, c("utm_x_scale", "utm_y_scale"), drop = FALSE],
      tid_values = tid_values
    )
    sm_pop_proj <- .normalize_smoother_output(
      parse_smoothers(NULL, data_utm, newdata = projection_data_use),
      n_rows = nrow(key) * n_t
    )
    fx_pop_proj <- .normalize_fixed_output(
      parse_fixed_covariates(NULL, data_utm, newdata = projection_data_use),
      n_rows = nrow(key) * n_t
    )
    rf_pop_proj <- .normalize_random_factor_output(
      parse_random_factors(NULL, data_utm, newdata = projection_data_use),
      n_rows = nrow(key) * n_t
    )
  }
  
  # ---- assemble data list for TMB ----
  data <- list(
    n_i = n_i,
    n_t = n_t,
    n_f = n_f,
    n_g = nrow(key),
    n_m = month_prep$n_m,
    
    b_i = data_utm$cpue,
    t_i = t_i,
    f_i = f_i,
    month_i = month_prep$month_i,
    has_tf = has_tf,
    area_g = key$area_km2_scaled,
    use_month_fe = month_prep$use_month_fe,
    
    A_is   = A_is,
    A_gs   = A_gs,
    Ais_ij = Ais_ij,
    Ais_x  = Ais_x,

    has_smooths_catch = as.integer(isTRUE(sm_catch$has_smooths)),
    Xs_catch = sm_catch$Xs,
    Zs_catch = sm_catch$Zs,
    b_smooth_start_catch = as.integer(sm_catch$b_smooth_start),
    has_fixed_catch = as.integer(isTRUE(fx_catch$has_fixed)),
    Xf_catch = fx_catch$X,
    has_random_factors_catch = as.integer(isTRUE(rf_catch$has_random)),
    Zf_catch = rf_catch$Z,
    b_factor_start_catch = as.integer(rf_catch$rand_start),

    has_smooths_pop = as.integer(isTRUE(sm_pop_obs$has_smooths)),
    Xs_pop_i = sm_pop_obs$Xs,
    Zs_pop_i = sm_pop_obs$Zs,
    Xs_pop_g = sm_pop_proj$Xs,
    Zs_pop_g = sm_pop_proj$Zs,
    b_smooth_start_pop = as.integer(sm_pop_obs$b_smooth_start),
    has_fixed_pop = as.integer(isTRUE(fx_pop_obs$has_fixed)),
    Xf_pop_i = fx_pop_obs$X,
    Xf_pop_g = fx_pop_proj$X,
    has_random_factors_pop = as.integer(isTRUE(rf_pop_obs$has_random)),
    Zf_pop_i = rf_pop_obs$Z,
    Zf_pop_g = rf_pop_proj$Z,
    b_factor_start_pop = as.integer(rf_pop_obs$rand_start),
    
    # PC priors
    matern_range   = diff(range(mesh$loc[, 1])) / 5,
    range_prob     = 0.5,
    matern_sigma_0 = 1,
    matern_sigma_t = 1,
    matern_sigma_fleet = 1,
    sigma_prob     = 0.05
  )
  
  data$spde <- .prep_anisotropy(mesh = mesh, spde = spde)
  
  list(
    data = data,
    key = key,
    scales = list(area_scale = key_out$area_scale_val),
    projection_data = projection_data_use,
    smooth_basis = list(
      catchability = sm_catch$basis_out,
      population = sm_pop_obs$basis_out
    ),
    covariate_info = list(
      catchability = list(
        smooth_labels = sm_catch$labels,
        K_smooth = ncol(sm_catch$Xs),
        n_smooth = length(sm_catch$Zs),
        sm_dims = sm_catch$sm_dims,
        fixed_terms = fx_catch$basis_out,
        K_fixed = ncol(fx_catch$X),
        random_factor_terms = rf_catch$basis_out,
        n_factor_random = length(rf_catch$Z),
        factor_rand_dims = rf_catch$rand_dims
      ),
      population = list(
        smooth_labels = sm_pop_obs$labels,
        K_smooth = ncol(sm_pop_obs$Xs),
        n_smooth = length(sm_pop_obs$Zs),
        sm_dims = sm_pop_obs$sm_dims,
        fixed_terms = fx_pop_obs$basis_out,
        K_fixed = ncol(fx_pop_obs$X),
        random_factor_terms = rf_pop_obs$basis_out,
        n_factor_random = length(rf_pop_obs$Z),
        factor_rand_dims = rf_pop_obs$rand_dims
      )
    ),
    time = list(
      index = index,
      month_diffs = month_diffs,
      month_col = month_prep$month_col,
      month_values = month_prep$month_values,
      use_month_fe = month_prep$use_month_fe
    )
  )
}

.prep_month_index <- function(
    data_utm,
    index = c("monthly", "yearly"),
    month_col = NULL,
    month_diffs = c("on", "off")
) {
  index <- match.arg(index)
  month_diffs <- match.arg(month_diffs)
  n_i <- nrow(data_utm)
  
  if (identical(index, "monthly") || identical(month_diffs, "off")) {
    return(list(
      n_m = 1L,
      month_i = rep.int(0L, n_i),
      month_values = integer(0),
      month_col = NULL,
      use_month_fe = 0L
    ))
  }
  
  if (is.null(month_col)) {
    candidates <- c("month", "Month")
    month_col <- candidates[candidates %in% names(data_utm)][1]
  }
  
  if (is.na(month_col) || is.null(month_col) || !month_col %in% names(data_utm)) {
    stop(
      paste0(
        "`index = \"yearly\"` with `month_diffs = \"on\"` requires a calendar month column. ",
        "Add `month` or `Month`, or supply `month_col`."
      ),
      call. = FALSE
    )
  }
  
  month_num <- suppressWarnings(as.numeric(as.character(data_utm[[month_col]])))
  
  if (anyNA(month_num) || any(!is.finite(month_num))) {
    stop("Calendar month values must be numeric and must not contain NA.", call. = FALSE)
  }
  if (any(abs(month_num - round(month_num)) > 1e-9)) {
    stop("Calendar month values must be integers in 1:12.", call. = FALSE)
  }
  
  month_int <- as.integer(round(month_num))
  if (any(month_int < 1L | month_int > 12L)) {
    stop("Calendar month values must be in 1:12.", call. = FALSE)
  }
  # Use only months observed in this dataset (e.g. left 170 may have 1,2,4..12; right 170 may have 4..11).
  # Sum-to-zero month FE is over these levels; each dataset gets its own n_m and month_i.
  month_values <- sort(unique(month_int))
  month_i <- match(month_int, month_values) - 1L
  
  list(
    n_m = length(month_values),
    month_i = as.integer(month_i),
    month_values = month_values,
    month_col = month_col,
    use_month_fe = as.integer(length(month_values) > 1L)
  )
}
