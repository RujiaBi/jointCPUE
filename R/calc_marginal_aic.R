#' Calculate marginal AIC
#'
#' Calculates an approximate marginal AIC by removing user-defined prior
#' penalties stored in `nll_prior` from the optimized objective value.
#'
#' @param object An object of class `jointCPUE` returned by
#'   [jointCPUE::jointCPUE()].
#'
#' @return A one-row data.frame with columns `logLik_marginal`, `k`,
#'   and `AIC_marginal`.
#' @export
calc_marginal_aic <- function(object) {
  if (!inherits(object, "jointCPUE")) {
    stop("`object` must be an `jointCPUE` fit from jointCPUE().", call. = FALSE)
  }

  par_best <- object$obj$env$last.par.best
  if (is.null(par_best)) {
    stop("Could not find `last.par.best` in fitted object.", call. = FALSE)
  }
  rep_obj <- object$obj$report(par_best)

  nll_total <- object$opt$objective
  nll_prior <- rep_obj$nll_prior
  nll_marginal <- nll_total - nll_prior

  k <- length(object$opt$par)

  data.frame(
    logLik_marginal = -nll_marginal,
    k = k,
    AIC_marginal = 2 * nll_marginal + 2 * k
  )
}
