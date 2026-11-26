\name{load_dataset}
\alias{load_dataset}
\title{Load a bundled example dataset}
\description{
Returns a numeric vector of sample men's 100 m times (seconds) shipped
with the package, for testing and demonstrating the algorithms.
}
\usage{
load_dataset(name = NULL)
}
\arguments{
  \item{name}{\code{"athletics20"} (20 times) or
    \code{"worldathletics-mens-100m-2379"} (2,379 times).
    If \code{NULL}, the known names are returned.}
}
\details{
Source: World Athletics 100m all-time lists, retrieved January 2024.

Sample Data Attribution: The sample 100m timing datasets (top 20 performances and 2379 all time bests) contain historical athletic performance records compiled from publicly available statistics published by World Athletics (https://worldathletics.org/). These datasets are provided strictly as static benchmark/example inputs for testing and demonstrating the algorithms in this package. This project is not affiliated with, endorsed by, or sponsored by World Athletics.

Each file keeps only the time column. Athlete names and other fields from the original lists are not included. The 20-row file is the prefix of
the 2,379-row list.
}
\value{
A numeric vector, or a character vector of names if \code{name} is
\code{NULL}.
}
\examples{
load_dataset()
x <- load_dataset("athletics20")
length(x)
res <- cluster(x, 4L)
res$bounds
res$sse
}
