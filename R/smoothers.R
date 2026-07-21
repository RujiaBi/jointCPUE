# Smoother parsing helpers adapted from the intCPUE workflow.
# They convert mgcv::s() terms into fixed and penalized design matrices.

.expand_back <- function(M_cc, keep, n_all) {
  if (is.null(M_cc)) return(NULL)
  if (is.vector(M_cc)) M_cc <- matrix(M_cc, ncol = 1L)
  M_all <- matrix(0, nrow = n_all, ncol = ncol(M_cc))
  if (nrow(M_cc) > 0L && any(keep)) M_all[keep, ] <- M_cc
  M_all
}

.complete_rows_for_sm <- function(sm, data) {
  needed <- sm$term
  if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L &&
      nzchar(sm$by) && sm$by %in% names(data)) {
    needed <- unique(c(needed, sm$by))
  }
  keep <- rep(TRUE, nrow(data))
  for (nm in needed) {
    if (!nm %in% names(data)) {
      stop("Smooth term requires variable `", nm, "` but it is not in data.", call. = FALSE)
    }
    keep <- keep & !is.na(data[[nm]])
  }
  keep
}

rm_wsp <- function(x) {
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}

all_terms <- function(x) {
  if (is.null(x) || !length(x)) {
    return(character(0))
  }
  if (!inherits(x, "formula")) {
    x <- stats::as.formula(x)
  }
  rm_wsp(attr(stats::terms(x), "term.labels"))
}

get_smooth_terms <- function(terms) {
  grep("s\\(", terms)
}

s2rPred <- function(sm, re, data) {
  needed <- sm$term
  if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L &&
      nzchar(sm$by) && sm$by %in% colnames(data)) {
    needed <- unique(c(needed, sm$by))
  }
  if (!all(needed %in% colnames(data))) {
    miss <- needed[!needed %in% colnames(data)]
    stop(
      "A smoother term is missing from `newdata`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  X <- mgcv::PredictMat(sm, data)
  if (!is.null(re$trans.U)) {
    X <- X %*% re$trans.U
  }
  X <- t(t(X) * re$trans.D)
  X[, re$rind] <- X[, re$pen.ind != 0]
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]

  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for (i in seq_along(re$rand)) {
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(r$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(r$rand) <- names(re$rand)
  r
}

parse_smoothers <- function(formula, data, knots = NULL, newdata = NULL, basis_prev = NULL) {
  terms <- all_terms(formula)
  smooth_i <- get_smooth_terms(terms)

  basis <- list()
  basis_out <- list()
  Zs <- list()
  Xs <- list()
  labels <- list()
  classes <- list()

  eval_data <- if (is.null(newdata)) data else newdata
  n_all <- nrow(eval_data)

  if (length(smooth_i) > 0L) {
    smterms <- terms[smooth_i]
    ns <- 0L
    ns_Xf <- 0L

    for (i in seq_along(smterms)) {
      if (grepl('bs\\=\\"re', smterms[i])) {
        stop("bs = 're' is not currently supported for smooths.", call. = FALSE)
      }
      if (grepl("fx\\=T", smterms[i]) || grepl("fx\\=TRUE", smterms[i])) {
        stop("fx = TRUE is not currently supported for smooths.", call. = FALSE)
      }

      expr <- str2expression(smterms[i])[[1]]
      eval_env <- new.env(parent = baseenv())
      eval_env$s <- mgcv::s
      obj <- eval(expr, envir = eval_env)

      labels[[i]] <- obj$label
      classes[[i]] <- attr(obj, "class")
      keep_i <- .complete_rows_for_sm(obj, eval_data)

      if (is.null(newdata)) {
        data_cc <- eval_data[keep_i, , drop = FALSE]
        if (nrow(data_cc) < 5L) {
          stop(
            "Too few non-missing rows (", nrow(data_cc),
            ") to build smooth term: ", labels[[i]],
            call. = FALSE
          )
        }
        basis_out[[i]] <- basis[[i]] <- mgcv::smoothCon(
          object = obj,
          data = data_cc,
          knots = knots,
          absorb.cons = TRUE,
          modCon = 3,
          diagonal.penalty = FALSE
        )
      } else {
        if (is.null(basis_prev) || length(basis_prev) < length(smterms)) {
          stop("`basis_prev` must have at least one entry per smooth term.", call. = FALSE)
        }
        basis[[i]] <- basis_prev[[i]]
      }

      for (j in seq_along(basis[[i]])) {
        ns_Xf <- ns_Xf + 1L

        if (is.null(newdata)) {
          data_cc <- eval_data[keep_i, , drop = FALSE]
          rasm_cc <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
          Xf <- .expand_back(rasm_cc$Xf, keep_i, n_all)
          rand_list <- lapply(rasm_cc$rand, .expand_back, keep = keep_i, n_all = n_all)
        } else {
          new_cc <- eval_data[keep_i, , drop = FALSE]
          rasm0 <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
          rasm_cc <- s2rPred(basis[[i]][[j]], rasm0, new_cc)
          Xf <- .expand_back(rasm_cc$Xf, keep_i, n_all)
          rand_list <- lapply(rasm_cc$rand, .expand_back, keep = keep_i, n_all = n_all)
        }

        for (k in seq_along(rand_list)) {
          if (is.null(rand_list[[k]])) next
          ns <- ns + 1L
          Zs[[ns]] <- Matrix::Matrix(rand_list[[k]], sparse = TRUE)
        }
        Xs[[ns_Xf]] <- Xf
      }
    }

    sm_dims <- as.integer(vapply(Zs, ncol, integer(1)))
    Xs <- if (length(Xs)) do.call(cbind, Xs) else matrix(0, nrow = n_all, ncol = 0L)
    b_smooth_start <- if (length(sm_dims)) c(0L, cumsum(sm_dims)[-length(sm_dims)]) else integer(0)
    has_smooths <- TRUE
  } else {
    has_smooths <- FALSE
    sm_dims <- integer(0)
    b_smooth_start <- integer(0)
    Xs <- matrix(0, nrow = n_all, ncol = 0L)
  }

  list(
    Xs = Xs,
    Zs = Zs,
    has_smooths = has_smooths,
    labels = labels,
    classes = classes,
    basis_out = basis_out,
    sm_dims = sm_dims,
    b_smooth_start = b_smooth_start
  )
}

.normalize_smoother_output <- function(sm, n_rows) {
  if (!isTRUE(sm$has_smooths)) {
    return(list(
      has_smooths = FALSE,
      Xs = matrix(0, nrow = n_rows, ncol = 0L),
      Zs = list(),
      basis_out = list(),
      labels = list(),
      classes = list(),
      sm_dims = integer(0),
      b_smooth_start = integer(0)
    ))
  }

  if (!identical(nrow(sm$Xs), n_rows)) {
    stop("parse_smoothers() returned `Xs` with unexpected row count.", call. = FALSE)
  }
  if (length(sm$Zs)) {
    for (k in seq_along(sm$Zs)) {
      if (!identical(nrow(sm$Zs[[k]]), n_rows)) {
        stop("parse_smoothers() returned a `Zs` matrix with unexpected row count.", call. = FALSE)
      }
      if (!inherits(sm$Zs[[k]], "sparseMatrix")) {
        sm$Zs[[k]] <- Matrix::Matrix(sm$Zs[[k]], sparse = TRUE)
      }
    }
  }

  list(
    has_smooths = TRUE,
    Xs = sm$Xs,
    Zs = sm$Zs,
    basis_out = sm$basis_out,
    labels = sm$labels,
    classes = sm$classes,
    sm_dims = as.integer(sm$sm_dims),
    b_smooth_start = as.integer(sm$b_smooth_start)
  )
}

.smooth_vars_from_basis <- function(basis_out) {
  needed <- character(0)
  if (!length(basis_out)) return(needed)

  keep_var_names <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & nzchar(x) & x != "NA"]
  }

  for (basis_i in basis_out) {
    if (!length(basis_i)) next
    for (sm in basis_i) {
      if (length(sm$term)) needed <- c(needed, keep_var_names(sm$term))
      if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L) {
        needed <- c(needed, keep_var_names(sm$by))
      }
    }
  }

  unique(keep_var_names(needed))
}
