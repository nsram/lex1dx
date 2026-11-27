\name{cluster}
\alias{cluster}
\title{Cluster 1-D data with lex1dx}
\description{
Partitions a numeric vector into \code{k} contiguous bins of the sorted
values, using the lex1dx pipeline, optional significant-digit grouping,
or an exact solve on the atoms.
}
\usage{
cluster(x, k, mode = "tuned", sig = NULL, M = NULL, tolerance = 1e-11,
        weights = NULL, solver = "wilber")
}
\arguments{
  \item{x}{A numeric vector. Empty or non-finite values are rejected.}
  \item{k}{A whole number, the requested number of bins.}
  \item{mode}{\code{"tuned"} (the default pipeline), \code{"sig5"}
    (an alias for \code{sig = 5}), or \code{"exact"} (no Collapse;
    the chosen solver runs on the atoms).}
  \item{sig}{An integer in 0..9, or \code{NULL}. 0 is full precision.
    1..9 groups to that many significant digits. If that grouping
    yields fewer than \code{k} groups, precision is raised (up to 15,
    then raw atoms) so \code{k} is kept when the raw data has at least
    \code{k} distinct values. \code{sig} in the result is the value
    actually used (0 if the fallback reached raw atoms).
    \code{NULL} follows \code{mode}. Values of 10 or more are
    rejected; use \code{sig = 0} for full precision.}
  \item{M}{Collapse target. \code{NULL} (the default) uses \code{default_M(k)}.
    Ignored when \code{mode} is \code{"exact"}.}
  \item{tolerance}{Refine move tolerance, a number in (0, 1), passed to the
    C++ core.}
  \item{weights}{Optional non-negative finite weights, the same length as
    \code{x}. \code{NULL} is the unweighted solver. The result SSE and
    centers are then weighted; \code{sizes} remain element counts.}
  \item{solver}{\code{"wilber"} (the default, modified Grønlund
    wilber-interp DP) or \code{"ckmeans"} (Fisher / ckmeans exact DP).
    Both solvers are exact and, when Wilber returns, give the same
    optimum (up to floating-point ties). On some inputs Wilber's
    \eqn{\lambda}{lambda} search can fail (an empty interval). The cause
    is collinear SSE*(k), e.g. regular grids or equally spaced values
    with equal counts, not ties. lex1dx then uses the Fisher / ckmeans
    solver, which is also exact. If that fallback's back-pointer table
    would exceed 2 GiB (possible in mode \code{"exact"} when D is large),
    lex1dx does not allocate the table and raises an error suggesting
    mode \code{"tuned"} or a smaller k. \code{solver_used} in the result
    is \code{"wilber"} or \code{"ckmeans"} and records which DP actually
    ran (when the DP is skipped, it mirrors the requested solver).
    Wilber is faster as k grows, because wilber-interp is O(k N) in
    style while ckmeans is a divide-and-conquer O(k N log N) DP.
    ckmeans is the original LEX exact DP and serves as a cross-check.}
}
\value{
A list with \code{bounds}, \code{sizes}, \code{centers}, \code{upper},
\code{sse}, \code{D}, \code{Dprime}, \code{sig}, \code{sweeps},
\code{moves}, \code{presorted}, and \code{solver_used}.
}
\seealso{
  \code{\link{default_M}}, \code{\link{load_dataset}},
  \code{\link{histogram}}, \code{\link{bin_table}},
  and the package vignette \code{vignette("lex1dx")}.
}
\examples{
x <- load_dataset("athletics20")
res <- cluster(x, 4L)
res$bounds
res$sse
res_sig <- cluster(x, 4L, sig = 5L)
res_sig$sig
vals <- unique(x)
wts <- as.numeric(table(x))
res_w <- cluster(vals, 4L, weights = wts)
res_w$sse
res_ex <- cluster(x, 4L, mode = "exact", solver = "ckmeans")
res_ex$sse
}
