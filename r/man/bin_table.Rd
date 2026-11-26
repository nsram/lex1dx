\name{bin_table}
\alias{bin_table}
\title{Table of lex1dx bins}
\description{
Returns one row per fitted bin: the minimum and maximum raw values,
the count or weight mass, the percentage of the total, and the
(weighted) mean.
}
\usage{
bin_table(result, x, weights = NULL)
}
\arguments{
  \item{result}{A result from \code{\link{cluster}}.}
  \item{x}{The numeric vector that was clustered.}
  \item{weights}{Optional non-negative weights, the same length as
    \code{x}.}
}
\value{
A data frame with columns \code{lo}, \code{hi}, \code{n}, \code{pct},
and \code{mean}.
}
\seealso{
  \code{\link{cluster}}, \code{\link{histogram}}.
}
\examples{
x <- load_dataset("athletics20")
r <- cluster(x, 5L)
bin_table(r, x)
}
