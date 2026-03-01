# Internal: build the TMB map for the jointCPUE model
.make_map_jointCPUE <- function(parameters, n_f,
                              pop_spatial = c("on", "off"),
                              pop_spatiotemporal = c("on", "off"),
                              q_diffs_system = c("on", "off"),
                              q_diffs_time = c("on", "off"),
                              q_diffs_spatial = c("on", "off"),
                              has_tf = NULL) {
  pop_spatial <- match.arg(pop_spatial)
  pop_spatiotemporal <- match.arg(pop_spatiotemporal)
  q_diffs_system <- match.arg(q_diffs_system)
  q_diffs_time <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)

  map <- list()

  if (pop_spatial == "off") {
    map$omega_s_1 <- factor(rep(NA, length(parameters$omega_s_1)))
    map$ln_sigma_0_1 <- factor(NA)
  }

  if (pop_spatiotemporal == "off") {
    map$epsilon_st_1 <- .map_matrix_NA(parameters$epsilon_st_1)
    map$ln_sigma_t_1 <- factor(NA)
  }

  if (q_diffs_system == "off" || n_f <= 1L) {
    map$fleet_f <- factor(rep(NA, length(parameters$fleet_f)))
    map$fleet_ln_std_dev <- factor(NA)
  }

  if (q_diffs_time == "off" || n_f <= 1L) {
    map$fleet_t <- .map_matrix_NA(parameters$fleet_t)
    map$fleet_t_ln_std_dev <- factor(NA)
  } else {
    if (!is.null(has_tf)) {
      keep <- has_tf > 0
      if (!all(dim(parameters$fleet_t) == dim(keep))) {
        stop("`has_tf` must have the same dimensions as `parameters$fleet_t`.", call. = FALSE)
      }
      map$fleet_t <- .map_matrix_partial_fix(parameters$fleet_t, keep)
      if (!any(keep)) {
        map$fleet_t_ln_std_dev <- factor(NA)
      }
    }
  }

  if (q_diffs_spatial == "off" || n_f <= 1L) {
    map$fleet_s <- factor(rep(NA, length(parameters$fleet_s)))
    map$ln_sigma_fleet <- factor(NA)
  }

  use_any_spde <- (
    pop_spatial == "on" ||
    pop_spatiotemporal == "on" ||
    (q_diffs_spatial == "on" && n_f > 1L)
  )

  if (!use_any_spde) {
    map$ln_range_1 <- factor(NA)
    map$ln_H_input <- factor(rep(NA, length(parameters$ln_H_input)))
  }

  map
}
