default_M <- function(k) {
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 1 || k != round(k) || k > .Machine$integer.max)
    stop("k must be a whole number >= 1")
  .default_M(as.integer(k))
}

cluster <- function(x, k, mode = "tuned", M = NULL, tolerance = 1e-11) {
  if (!is.numeric(x) || (!is.null(dim(x)) && length(dim(x)) > 1L))
    stop("x must be a numeric vector")
  if (length(x) == 0L) stop("x is empty")
  if (anyNA(x)) stop("x contains NA or NaN")
  if (any(is.infinite(x))) stop("x contains Inf or -Inf")
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 1 || k != round(k) || k > .Machine$integer.max)
    stop("k must be a whole number >= 1")
  k <- as.integer(k)
  if (!is.character(mode) || length(mode) != 1L || !mode %in% c("tuned", "sig5"))
    stop("mode must be 'tuned' or 'sig5'")
  if (is.null(M)) {
    coarse <- 0
  } else {
    if (length(M) != 1L || !is.numeric(M) || !is.finite(M) || M < 0 || M != round(M))
      stop("M must be a whole number >= 0")
    coarse <- as.double(M)
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) || tolerance < 0)
    stop("tolerance must be a finite number >= 0")
  .cluster(as.double(x), k, mode, coarse, as.double(tolerance))
}
