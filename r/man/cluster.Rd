\name{cluster}
\alias{cluster}
\title{Cluster 1-D data with lex1dx}
\description{
Partitions a numeric vector into \code{k} contiguous bins of the sorted
values, using the lex1dx tuned pipeline or optional 5-significant-digit grouping.
}
\usage{
cluster(x, k, mode = "tuned", M = NULL, tolerance = 1e-11)
}
\arguments{
  \item{x}{A numeric vector. Empty or non-finite values are rejected.}
  \item{k}{A whole number, the requested number of bins.}
  \item{mode}{\code{"tuned"} (full-precision atoms) or \code{"sig5"}.}
  \item{M}{Collapse target. \code{NULL} (the default) uses \code{default_M(k)}.}
  \item{tolerance}{Refine move tolerance, passed to the C++ core.}
}
\value{
A list with \code{bounds}, \code{sizes}, \code{centers}, \code{upper},
\code{sse}, \code{D}, \code{Dprime}, \code{sweeps}, \code{moves}, and
\code{presorted}.
}
