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
    month_diffs = c("on", "off")
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
  key_out <- .prep_key_area(data_utm, mesh, area_scale = area_scale)
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
  
  t_i <- t_chk$x
  f_i <- f_chk$x
  
  n_t <- t_chk$n
  n_f <- f_chk$n

  has_tf <- matrix(0L, nrow = n_t, ncol = max(0L, n_f - 1L))
  if (n_f > 1L) {
    ii <- which(f_i > 0L)
    if (length(ii) > 0L) {
      tt <- t_i[ii] + 1L
      ff <- f_i[ii]
      has_tf[cbind(tt, ff)] <- 1L
    }
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
