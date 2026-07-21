.extract_formula_terms <- function(formula) {
  terms <- all_terms(formula)
  if (!length(terms)) {
    return(list(all = character(0), smooth = character(0), random = character(0), fixed = character(0)))
  }
  smooth_i <- get_smooth_terms(terms)
  random_i <- grep("^f\\(", terms)
  non_fixed_i <- unique(c(smooth_i, random_i))
  list(
    all = terms,
    smooth = terms[smooth_i],
    random = terms[random_i],
    fixed = if (length(non_fixed_i)) terms[-non_fixed_i] else terms
  )
}

.parse_simple_term_var <- function(term) {
  if (!grepl("^[[:alnum:]_.]+$", term)) {
    stop(
      "Only simple fixed covariate terms like `season` are supported outside `s()`. ",
      "Unsupported term: ", term,
      call. = FALSE
    )
  }
  term
}

.parse_random_factor_var <- function(term) {
  m <- regexec("^f\\(([[:alnum:]_.]+)\\)$", term)
  reg <- regmatches(term, m)[[1]]
  if (length(reg) != 2L) {
    stop(
      "Random factor terms must use simple syntax like `f(season)`. ",
      "Unsupported term: ", term,
      call. = FALSE
    )
  }
  reg[2]
}

.normalize_factor_levels <- function(x, term) {
  x_chr <- as.character(x)
  x_chr <- x_chr[!is.na(x_chr)]
  if (!length(x_chr)) {
    stop("Factor term `", term, "` has no non-missing values.", call. = FALSE)
  }
  if (is.factor(x)) {
    lvls <- levels(x)
    lvls[lvls %in% unique(x_chr)]
  } else {
    sort(unique(x_chr))
  }
}

.build_sum_to_zero_factor_matrix <- function(values, levels, term) {
  fac <- factor(as.character(values), levels = levels)
  n_level <- length(levels)
  if (n_level < 2L) {
    return(matrix(0, nrow = length(values), ncol = 0L))
  }
  mm <- matrix(0, nrow = length(values), ncol = n_level - 1L)
  idx <- as.integer(fac)
  for (j in seq_len(n_level - 1L)) {
    mm[idx == j, j] <- 1
  }
  mm[idx == n_level, ] <- -1
  colnames(mm) <- paste0(term, "::", levels[-n_level])
  mm
}

.build_linear_matrix <- function(values, term) {
  x <- suppressWarnings(as.numeric(values))
  if (any(!is.finite(x))) {
    stop("Numeric covariate `", term, "` must be finite and must not contain NA.", call. = FALSE)
  }
  out <- matrix(x, ncol = 1L)
  colnames(out) <- term
  out
}

.expand_back_sparse <- function(M_cc, keep, n_all) {
  if (is.null(M_cc)) return(NULL)
  if (!inherits(M_cc, "sparseMatrix")) {
    M_cc <- Matrix::Matrix(M_cc, sparse = TRUE)
  }
  M_all <- Matrix::Matrix(0, nrow = n_all, ncol = ncol(M_cc), sparse = TRUE)
  if (nrow(M_cc) > 0L && any(keep)) {
    M_all[keep, ] <- M_cc
  }
  M_all
}

.build_random_factor_matrix <- function(values, levels, term) {
  fac <- factor(as.character(values), levels = levels)
  mm <- Matrix::sparse.model.matrix(~ fac - 1)
  colnames(mm) <- paste0(term, "::", levels)
  mm
}

parse_fixed_covariates <- function(formula, data, newdata = NULL, basis_prev = NULL) {
  parts <- .extract_formula_terms(formula)
  eval_data <- if (is.null(newdata)) data else newdata
  n_all <- nrow(eval_data)

  if (!length(parts$fixed)) {
    return(list(
      X = matrix(0, nrow = n_all, ncol = 0L),
      has_fixed = FALSE,
      basis_out = list()
    ))
  }

  X_list <- list()
  basis_out <- list()

  for (i in seq_along(parts$fixed)) {
    term <- parts$fixed[i]
    var <- .parse_simple_term_var(term)
    if (!var %in% names(eval_data)) {
      stop("Fixed covariate term `", term, "` requires variable `", var, "` in data.", call. = FALSE)
    }

    if (is.null(newdata)) {
      values <- eval_data[[var]]
      kind <- if (is.numeric(values) || is.integer(values)) "numeric" else "factor"
      if (identical(kind, "numeric")) {
        basis_i <- list(term = term, var = var, kind = "numeric")
        X_i <- .build_linear_matrix(values, term)
      } else {
        levels_i <- .normalize_factor_levels(values, term)
        basis_i <- list(term = term, var = var, kind = "factor", levels = levels_i)
        X_i <- .build_sum_to_zero_factor_matrix(values, levels_i, term)
      }
    } else {
      if (is.null(basis_prev) || length(basis_prev) < i) {
        stop("`basis_prev` is missing fixed covariate basis information for `", term, "`.", call. = FALSE)
      }
      basis_i <- basis_prev[[i]]
      values <- eval_data[[var]]
      if (identical(basis_i$kind, "numeric")) {
        X_i <- .build_linear_matrix(values, term)
      } else {
        vals_cc <- as.character(values[!is.na(values)])
        unseen <- setdiff(unique(vals_cc), basis_i$levels)
        if (length(unseen)) {
          stop(
            "Factor term `", term, "` in projection data contains unseen level(s): ",
            paste(unseen, collapse = ", "),
            call. = FALSE
          )
        }
        X_i <- .build_sum_to_zero_factor_matrix(values, basis_i$levels, term)
      }
    }

    X_list[[length(X_list) + 1L]] <- X_i
    basis_out[[i]] <- basis_i
  }

  X <- if (length(X_list)) do.call(cbind, X_list) else matrix(0, nrow = n_all, ncol = 0L)
  list(
    X = X,
    has_fixed = ncol(X) > 0L,
    basis_out = basis_out
  )
}

.normalize_fixed_output <- function(fx, n_rows) {
  if (!identical(nrow(fx$X), n_rows)) {
    stop("parse_fixed_covariates() returned `X` with unexpected row count.", call. = FALSE)
  }
  list(
    X = fx$X,
    has_fixed = isTRUE(fx$has_fixed),
    basis_out = fx$basis_out
  )
}

.fixed_vars_from_basis <- function(basis_out) {
  if (!length(basis_out)) return(character(0))
  vars <- vapply(basis_out, `[[`, "", "var")
  unique(vars[!is.na(vars) & nzchar(vars)])
}

parse_random_factors <- function(formula, data, newdata = NULL, basis_prev = NULL) {
  parts <- .extract_formula_terms(formula)
  eval_data <- if (is.null(newdata)) data else newdata
  n_all <- nrow(eval_data)

  if (!length(parts$random)) {
    return(list(
      Z = list(),
      has_random = FALSE,
      basis_out = list(),
      rand_dims = integer(0),
      rand_start = integer(0)
    ))
  }

  Z_list <- list()
  basis_out <- list()

  for (i in seq_along(parts$random)) {
    term <- parts$random[i]
    var <- .parse_random_factor_var(term)
    if (!var %in% names(eval_data)) {
      stop("Random factor term `", term, "` requires variable `", var, "` in data.", call. = FALSE)
    }
    keep <- !is.na(eval_data[[var]])

    if (is.null(newdata)) {
      levels_i <- .normalize_factor_levels(eval_data[[var]], term)
      basis_i <- list(term = term, var = var, kind = "random_factor", levels = levels_i)
    } else {
      if (is.null(basis_prev) || length(basis_prev) < i) {
        stop("`basis_prev` is missing random factor basis information for `", term, "`.", call. = FALSE)
      }
      basis_i <- basis_prev[[i]]
      vals_cc <- as.character(eval_data[[var]][keep])
      unseen <- setdiff(unique(vals_cc), basis_i$levels)
      if (length(unseen)) {
        stop(
          "Random factor term `", term, "` in projection data contains unseen level(s): ",
          paste(unseen, collapse = ", "),
          call. = FALSE
        )
      }
    }

    Z_cc <- .build_random_factor_matrix(eval_data[[var]][keep], basis_i$levels, term)
    Z_list[[length(Z_list) + 1L]] <- .expand_back_sparse(Z_cc, keep, n_all)
    basis_out[[i]] <- basis_i
  }

  rand_dims <- as.integer(vapply(Z_list, ncol, integer(1)))
  rand_start <- if (length(rand_dims)) c(0L, cumsum(rand_dims)[-length(rand_dims)]) else integer(0)

  list(
    Z = Z_list,
    has_random = length(Z_list) > 0L,
    basis_out = basis_out,
    rand_dims = rand_dims,
    rand_start = rand_start
  )
}

.normalize_random_factor_output <- function(rf, n_rows) {
  if (length(rf$Z)) {
    for (k in seq_along(rf$Z)) {
      if (!identical(nrow(rf$Z[[k]]), n_rows)) {
        stop("parse_random_factors() returned a `Z` matrix with unexpected row count.", call. = FALSE)
      }
      if (!inherits(rf$Z[[k]], "sparseMatrix")) {
        rf$Z[[k]] <- Matrix::Matrix(rf$Z[[k]], sparse = TRUE)
      }
    }
  }

  list(
    Z = rf$Z,
    has_random = isTRUE(rf$has_random),
    basis_out = rf$basis_out,
    rand_dims = as.integer(rf$rand_dims),
    rand_start = as.integer(rf$rand_start)
  )
}

.random_factor_vars_from_basis <- function(basis_out) {
  if (!length(basis_out)) return(character(0))
  vars <- vapply(basis_out, `[[`, "", "var")
  unique(vars[!is.na(vars) & nzchar(vars)])
}

.formula_has_var <- function(formula, var) {
  !is.null(formula) && var %in% all.vars(formula)
}

.check_month_formula_conflict <- function(formula_catchability, formula_population, month_col, use_month_fe) {
  if (!isTRUE(use_month_fe == 1L) || is.null(month_col)) {
    return(invisible(TRUE))
  }
  if (.formula_has_var(formula_catchability, month_col) ||
      .formula_has_var(formula_population, month_col)) {
    stop(
      "`month_diffs = \"on\"` already adds `", month_col,
      "` as a sum-to-zero observation-only effect. Remove it from formula covariates.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.expand_projection_over_time <- function(projection_data, tid_values) {
  projection_data <- as.data.frame(projection_data)
  tid_values <- as.integer(tid_values)
  out <- projection_data[rep(seq_len(nrow(projection_data)), times = length(tid_values)), , drop = FALSE]
  out$tid <- rep(tid_values, each = nrow(projection_data))
  rownames(out) <- NULL
  out
}

.warn_projection_na <- function(projection_data, needed_vars) {
  if (!length(needed_vars)) {
    return(invisible(NULL))
  }
  na_counts <- vapply(
    needed_vars,
    function(nm) sum(is.na(projection_data[[nm]])),
    integer(1)
  )
  na_counts <- na_counts[na_counts > 0L]
  if (length(na_counts)) {
    warning(
      paste0(
        "`projection_data` contains NA in population covariates: ",
        paste(sprintf("%s (%d)", names(na_counts), as.integer(na_counts)), collapse = ", "),
        ". Affected smooth/factor rows contribute 0 only when their parser supports NA expansion."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.align_projection_data_to_key <- function(projection_data, key, needed_vars, tid_values) {
  projection_data <- as.data.frame(projection_data)
  tid_values <- as.integer(tid_values)
  coord_cols <- c("utm_x_scale", "utm_y_scale")
  has_tid <- "tid" %in% names(projection_data)
  req <- unique(c(coord_cols, if (has_tid) "tid", needed_vars))
  miss <- setdiff(req, names(projection_data))
  if (length(miss)) {
    stop("`projection_data` is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  if (has_tid) {
    if (anyNA(projection_data$tid)) {
      stop("`projection_data$tid` must not contain NA.", call. = FALSE)
    }
    tid_num <- suppressWarnings(as.numeric(as.character(projection_data$tid)))
    if (any(!is.finite(tid_num)) || any(abs(tid_num - round(tid_num)) > 1e-9)) {
      stop("`projection_data$tid` must be integer-valued and use the same coding as `data_utm$tid`.", call. = FALSE)
    }
    projection_data$tid <- as.integer(round(tid_num))
    bad_tid <- setdiff(unique(projection_data$tid), tid_values)
    if (length(bad_tid)) {
      stop(
        "`projection_data$tid` contains values not used by the model: ",
        paste(sort(bad_tid), collapse = ", "),
        call. = FALSE
      )
    }

    proj_id <- paste(projection_data$utm_x_scale, projection_data$utm_y_scale, projection_data$tid, sep = "\r")
    key_gt <- .expand_projection_over_time(key[, coord_cols, drop = FALSE], tid_values)
    key_id <- paste(key_gt$utm_x_scale, key_gt$utm_y_scale, key_gt$tid, sep = "\r")
    if (anyDuplicated(proj_id)) {
      stop("`projection_data` must contain at most one row per extrapolation cell-time combination.", call. = FALSE)
    }
    idx <- match(key_id, proj_id)
    if (anyNA(idx)) {
      stop(
        "`projection_data` must contain one row for every extrapolation grid cell-time combination.",
        call. = FALSE
      )
    }
    out <- projection_data[idx, req, drop = FALSE]
  } else {
    proj_id <- paste(projection_data$utm_x_scale, projection_data$utm_y_scale, sep = "\r")
    key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")
    if (anyDuplicated(proj_id)) {
      stop("`projection_data` must contain at most one row per extrapolation cell.", call. = FALSE)
    }
    idx <- match(key_id, proj_id)
    if (anyNA(idx)) {
      stop("`projection_data` must contain one row for every extrapolation grid cell.", call. = FALSE)
    }
    out <- .expand_projection_over_time(
      projection_data[idx, unique(c(coord_cols, needed_vars)), drop = FALSE],
      tid_values = tid_values
    )
  }

  rownames(out) <- NULL
  out
}

.default_population_projection_data <- function(data_utm, key, needed_vars, tid_values) {
  coord_cols <- c("utm_x_scale", "utm_y_scale")
  if (!length(needed_vars)) {
    return(.expand_projection_over_time(key[, coord_cols, drop = FALSE], tid_values))
  }

  miss <- setdiff(needed_vars, names(data_utm))
  if (length(miss)) {
    stop(
      "Population formula terms require covariates not found in `data_utm`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  cell_id <- paste(data_utm$utm_x_scale, data_utm$utm_y_scale, sep = "\r")
  key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")
  out_g <- key[, coord_cols, drop = FALSE]
  for (nm in needed_vars) {
    proto <- data_utm[[nm]]
    out_g[[nm]] <- proto[rep(NA_integer_, nrow(key))]
  }

  for (g in seq_len(nrow(key))) {
    rows <- which(cell_id == key_id[g])
    if (!length(rows)) {
      stop("Internal error: extrapolation cell not found in `data_utm`.", call. = FALSE)
    }
    for (nm in needed_vars) {
      vals <- data_utm[[nm]][rows]
      vals_non_na <- vals[!is.na(vals)]
      uniq <- unique(vals_non_na)
      if (length(uniq) > 1L) {
        stop(
          "Population covariate `", nm, "` varies within an extrapolation grid cell. ",
          "Please supply `projection_data` explicitly. Use a `tid` column there ",
          "if the covariate varies over time.",
          call. = FALSE
        )
      }
      if (length(uniq) == 1L) {
        out_g[[nm]][g] <- uniq
      }
    }
  }

  .expand_projection_over_time(out_g, tid_values = tid_values)
}
