\name{histogram}
\alias{histogram}
\title{Plot a lex1dx bin table}
\description{
Draws one horizontal bar per fitted bin, from the lowest bin at the
top. Each row shows the bin's minimum and maximum raw values, a bar
proportional to the count (or weight mass), and the count, percentage,
and mean. The header reports r-squared and MSE.
}
\usage{
histogram(result, x, file = NULL,
          weights = NULL, title = "", theme = "default")
}
\arguments{
  \item{result}{A result from \code{\link{cluster}}.}
  \item{x}{The numeric vector that was clustered.}
  \item{file}{A PNG path, or \code{NULL} (the default).
    \code{NULL} draws to the current device and writes no file.}
  \item{weights}{Optional non-negative weights, the same length as
    \code{x}. When given, bar length, \code{n}, and the mean use the
    weight mass.}
  \item{title}{An optional first header line. The default is empty.}
  \item{theme}{\code{"default"} is a blue background with yellow bars
    and text. \code{"light"} is a white background with dark grey bars
    and black text.}
}
\details{
r-squared is 1 - SSE/SST. MSE is SSE/N, or SSE divided by the total
weight when \code{weights} is given. Bin limits are shown to 5
significant digits. The plot is built from \code{\link{bin_table}}.
}
\value{
If \code{file} is not \code{NULL}, that path, invisibly. Otherwise
\code{NULL}, invisibly.
}
\seealso{
  \code{\link{cluster}}, \code{\link{bin_table}}.
}
\examples{
x <- load_dataset("athletics20")
r <- cluster(x, 5L)
f <- tempfile(fileext = ".png")
histogram(r, x, file = f, title = "athletics20")
}
