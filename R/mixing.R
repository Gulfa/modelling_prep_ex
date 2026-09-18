#' Combine an age mixing matrix and a region mixing matrix
#'
#' The backend model works on `n_age * n_region` groups, numbered with age
#' varying fastest: group `(r - 1) * n_age + a` is age `a` in region `r`. The
#' contact rate between age `a` in region `r` and age `b` in region `s` is
#'
#' `C[(a, r), (b, s)] = C_age[a, b] * M_region[r, s]`
#'
#' i.e. everyone has the age-mixing pattern `C_age`, and `M_region` decides in
#' which regions those contacts take place.
#'
#' @param C_age `n_age x n_age` matrix; `C_age[a, b]` is the number of contacts
#'   per day a person of age `a` has with people of age `b`.
#' @param M_region `n_region x n_region` matrix; `M_region[r, s]` is the
#'   fraction of a region-`r` resident's contacts made with residents of region
#'   `s`. Rows must sum to 1. The identity matrix gives isolated regions.
#'
#' @return The `n_age * n_region` square contact matrix.
#' @export
combine_mixing <- function(C_age, M_region) {
  kronecker(M_region, C_age)
}

check_C_age <- function(C_age, n_age, age_labels) {
  if (is.null(C_age)) {
    if (n_age > 1) {
      stop("'C_age' is required when there is more than one age group")
    }
    return(matrix(1, 1, 1, dimnames = list(age_labels, age_labels)))
  }
  C_age <- as.matrix(C_age)
  if (!identical(dim(C_age), c(n_age, n_age))) {
    stop(sprintf("'C_age' must be a %d x %d matrix (one row/column per age group)",
                 n_age, n_age))
  }
  if (anyNA(C_age) || any(C_age < 0)) {
    stop("'C_age' must be non-negative")
  }
  dimnames(C_age) <- list(age_labels, age_labels)
  C_age
}

check_M_region <- function(M_region, n_region, region_labels) {
  if (is.null(M_region)) {
    if (n_region > 1) {
      stop("'M_region' is required when there is more than one region")
    }
    return(matrix(1, 1, 1, dimnames = list(region_labels, region_labels)))
  }
  M_region <- as.matrix(M_region)
  if (!identical(dim(M_region), c(n_region, n_region))) {
    stop(sprintf("'M_region' must be a %d x %d matrix (one row/column per region)",
                 n_region, n_region))
  }
  if (anyNA(M_region) || any(M_region < 0)) {
    stop("'M_region' must be non-negative")
  }
  if (any(abs(rowSums(M_region) - 1) > 1e-6)) {
    stop("Rows of 'M_region' must sum to 1 (fraction of contacts made in each region)")
  }
  dimnames(M_region) <- list(region_labels, region_labels)
  M_region
}
