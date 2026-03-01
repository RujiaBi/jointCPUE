#' Get observation-level predicted values and residuals
#'
#' Extracts observation-level predicted values and residuals from a fitted
#' `jointCPUE` model. These can be used for residual diagnostics, including
#' spatial residual pattern plots.
#'
#' @param object An object of class `jointCPUE` returned by
#'   [jointCPUE::jointCPUE()].
#' @param data Optional original data.frame to bind to the returned fitted
#'   values. If supplied, it must have the same number of rows and order as the
#'   fitted observations.
#' @param drop_floor Logical. If `TRUE`, drops observations at or below a
#'   specified floor value. This is useful when a small constant (for example
#'   `1e-5`) has been added to zero CPUE values before fitting.
#' @param floor_value Optional numeric floor used when `drop_floor = TRUE`.
#'   If `NULL`, the minimum observed value in the fitted data is used.
#'
#' @return A data.frame with predicted values and residuals. If `data` is
#'   supplied, it is column-bound to the output.
#' @export
get_predicted <- function(object, data = NULL, drop_floor = FALSE, floor_value = NULL) {
  if (!inherits(object, "jointCPUE")) {
    stop("`object` must be an `jointCPUE` fit from jointCPUE().", call. = FALSE)
  }

  par_best <- object$obj$env$last.par.best
  if (is.null(par_best)) {
    stop("Could not find `last.par.best` in fitted object.", call. = FALSE)
  }
  rep_obj <- object$obj$report(par_best)
  req <- c("eta_hat_i", "mu_hat_i", "residual_raw_i", "residual_log_i")
  miss <- setdiff(req, names(rep_obj))
  if (length(miss)) {
    stop(
      "Missing fitted-value outputs in `obj$report()`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  out <- data.frame(
    obs = seq_along(rep_obj$mu_hat_i),
    observed = object$data_tmb$b_i,
    tid = object$data_tmb$t_i,
    fleetid = object$data_tmb$f_i,
    eta_hat = rep_obj$eta_hat_i,
    fitted = rep_obj$mu_hat_i,
    observed_log = log(object$data_tmb$b_i),
    predicted_log = rep_obj$eta_hat_i,
    residual_raw = rep_obj$residual_raw_i,
    residual_log = rep_obj$residual_log_i
  )

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame if supplied.", call. = FALSE)
    }
    if (nrow(data) != nrow(out)) {
      stop(
        "`data` must have the same number of rows as the fitted model. Got ",
        nrow(data), " rows but expected ", nrow(out), ".",
        call. = FALSE
      )
    }
    out <- cbind(data, out)
  }

  if (isTRUE(drop_floor)) {
    if (is.null(floor_value)) {
      floor_value <- min(out$observed, na.rm = TRUE)
    }
    if (!is.numeric(floor_value) || length(floor_value) != 1L || !is.finite(floor_value)) {
      stop("`floor_value` must be a single finite numeric value.", call. = FALSE)
    }
    out <- out[out$observed > floor_value, , drop = FALSE]
  }

  out
}

#' Get observation-level fitted values and residuals
#'
#' Backward-compatible wrapper for [get_predicted()].
#'
#' @inheritParams get_predicted
#' @return A data.frame with predicted values and residuals.
#' @export
get_fitted <- function(object, data = NULL, drop_floor = FALSE, floor_value = NULL) {
  get_predicted(
    object = object,
    data = data,
    drop_floor = drop_floor,
    floor_value = floor_value
  )
}
