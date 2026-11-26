default_M <- function(k) {
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 1 || k != round(k) || k > .Machine$integer.max)
    stop("k must be a whole number >= 1")
  .default_M(as.integer(k))
}

cluster <- function(x, k, mode = "tuned", sig = NULL, M = NULL, tolerance = 1e-11,
                    weights = NULL, solver = "wilber") {
  if (!is.numeric(x) || (!is.null(dim(x)) && length(dim(x)) > 1L))
    stop("x must be a numeric vector")
  if (length(x) == 0L) stop("x is empty")
  if (anyNA(x)) stop("x contains NA or NaN")
  if (any(is.infinite(x))) stop("x contains Inf or -Inf")
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 1 || k != round(k) || k > .Machine$integer.max)
    stop("k must be a whole number >= 1")
  k <- as.integer(k)
  if (!is.character(mode) || length(mode) != 1L || !mode %in% c("tuned", "sig5", "exact"))
    stop("mode must be 'tuned', 'sig5', or 'exact'")
  if (!is.character(solver) || length(solver) != 1L || !solver %in% c("wilber", "ckmeans"))
    stop("solver must be 'wilber' or 'ckmeans'")
  if (is.null(sig)) {
    sig_arg <- -1L
  } else {
    if (length(sig) != 1L || !is.numeric(sig) || !is.finite(sig) || sig != round(sig) || sig < 0 || sig > 9)
      stop("sig must be an integer in 0..9. Use sig=0 for full precision")
    sig_arg <- as.integer(sig)
  }
  if (is.null(M)) {
    coarse <- 0
  } else {
    if (length(M) != 1L || !is.numeric(M) || !is.finite(M) || M < 0 || M != round(M))
      stop("M must be a whole number >= 0")
    coarse <- as.double(M)
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) || tolerance < 0)
    stop("tolerance must be a finite number >= 0")
  w <- NULL
  if (!is.null(weights)) {
    if (!is.numeric(weights) || (!is.null(dim(weights)) && length(dim(weights)) > 1L))
      stop("weights must be a numeric vector")
    if (length(weights) != length(x))
      stop("weights size must match input")
    if (anyNA(weights))
      stop("weight is NA or NaN; weights must be finite and >= 0")
    if (any(is.infinite(weights)) || any(weights < 0))
      stop("weights must be finite and >= 0")
    if (!(sum(weights) > 0) || !is.finite(sum(weights)))
      stop("weights must sum to a positive finite total")
    w <- as.double(weights)
  }
  .cluster(as.double(x), k, mode, coarse, as.double(tolerance), sig_arg, w, solver)
}
