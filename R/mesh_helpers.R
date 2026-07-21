#' @keywords internal
.as_jointCPUEmesh <- function(mesh, loc_xy = NULL, xy_cols = NULL, recompute_A = c("auto", "always", "never")) {
  recompute_A <- match.arg(recompute_A)
  
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("Package 'fmesher' is required.", call. = FALSE)
  }
  
  # Case 1: already jointCPUEmesh
  if (inherits(mesh, "jointCPUEmesh")) {
    m <- mesh$mesh
    if (is.null(m)) stop("`mesh` is an jointCPUEmesh but `mesh$mesh` is NULL.", call. = FALSE)
    
    # Ensure spde exists
    spde <- mesh$spde
    if (is.null(spde)) spde <- fmesher::fm_fem(m)
    
    # Decide if A can be reused
    A <- mesh$A
    need_A <- FALSE
    
    if (is.null(loc_xy)) {
      # no loc requested -> keep existing A (if present)
      if (is.null(A)) need_A <- TRUE
    } else {
      # loc requested -> A must correspond to loc_xy
      if (is.null(A)) {
        need_A <- TRUE
      } else if (recompute_A == "always") {
        need_A <- TRUE
      } else if (recompute_A == "auto") {
        old_loc <- mesh$loc_xy
        same_locations <- !is.null(old_loc) &&
          identical(dim(as.matrix(old_loc)), dim(as.matrix(loc_xy))) &&
          identical(unname(as.matrix(old_loc)), unname(as.matrix(loc_xy)))
        if (!same_locations || !identical(nrow(A), nrow(loc_xy))) {
          need_A <- TRUE
        }
      }
    }
    
    if (need_A) {
      if (is.null(loc_xy)) stop("Need to compute A but `loc_xy` was not provided.", call. = FALSE)
      A <- fmesher::fm_basis(m, loc = loc_xy)
    }
    
    out <- mesh
    out$mesh <- m
    out$spde <- spde
    out$A <- A
    if (!is.null(loc_xy)) out$loc_xy <- loc_xy
    if (!is.null(xy_cols)) out$xy_cols <- xy_cols
    return(out)
  }
  
  # Case 2: bare fmesher mesh
  m <- mesh
  if (is.null(m)) stop("`mesh` is NULL.", call. = FALSE)
  
  if (is.null(loc_xy)) stop("When `mesh` is not an jointCPUEmesh, you must provide `loc_xy`.", call. = FALSE)
  
  spde <- fmesher::fm_fem(m)
  A <- fmesher::fm_basis(m, loc = loc_xy)
  
  structure(list(
    loc_xy = loc_xy,
    xy_cols = if (is.null(xy_cols)) c(NA_character_, NA_character_) else xy_cols,
    mesh = m,
    spde = spde,
    A = A,
    loc_centers = NULL
  ), class = "jointCPUEmesh")
}

.validate_fem_basis <- function(A, expected_rows, expected_cols, name, tolerance = 1e-6) {
  if (!inherits(A, "Matrix") && !is.matrix(A)) {
    stop("`", name, "` must be a matrix-like FEM basis object.", call. = FALSE)
  }
  if (!identical(nrow(A), as.integer(expected_rows)) ||
      !identical(ncol(A), as.integer(expected_cols))) {
    stop(
      "`", name, "` has dimensions ", nrow(A), " x ", ncol(A),
      "; expected ", expected_rows, " x ", expected_cols, ".",
      call. = FALSE
    )
  }

  values <- if (inherits(A, "sparseMatrix")) A@x else as.numeric(A)
  row_sums <- as.numeric(Matrix::rowSums(A))
  bad <- which(!is.finite(row_sums) | abs(row_sums - 1) > tolerance)

  if (any(!is.finite(values)) || length(bad)) {
    shown <- paste(utils::head(bad, 10L), collapse = ", ")
    suffix <- if (length(bad) > 10L) ", ..." else ""
    stop(
      "`", name, "` contains invalid FEM projection rows",
      if (length(bad)) paste0(" (rows: ", shown, suffix, ")") else "",
      ". Coordinates may lie outside the mesh triangulation.",
      call. = FALSE
    )
  }

  invisible(A)
}
