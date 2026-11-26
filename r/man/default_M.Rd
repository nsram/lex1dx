\name{default_M}
\alias{default_M}
\title{Default Collapse target}
\description{
Returns the default Collapse target for \code{k} bins: \code{max(4096, 32 * k)}.
}
\usage{
default_M(k)
}
\arguments{
  \item{k}{A whole number, the requested number of bins.}
}
\value{An integer.}
\examples{
default_M(64L)
default_M(256L)
}
