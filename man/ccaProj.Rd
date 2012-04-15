\name{ccaProj}
\alias{ccaProj}
\title{(Robust) CCA via projections through the data points}
\usage{
  ccaProj(x, y, k = 1,
    method = c("spearman", "kendall", "quadrant", "M", "pearson"),
    control = list(...), useL1Median = TRUE, fast = TRUE,
    initial = "random", nIterations = 10, seed = NULL, ...)
}
\arguments{
  \item{x,y}{each can be a numeric vector, matrix or data
  frame.}

  \item{k}{an integer giving the number of canonical
  variables to compute.}

  \item{method}{a character string specifying the
  correlation functional to maximize.  Possible values are
  \code{"spearman"} for the Spearman correlation,
  \code{"kendall"} for the Kendall correlation,
  \code{"quadrant"} for the quadrant correlation,
  \code{"M"} for the correlation based on a bivariate
  M-estimator of location and scatter with a Huber loss
  function, or \code{"pearson"} for the classical Pearson
  correlation (see \code{\link{corFunctions}}).}

  \item{control}{a list of additional arguments to be
  passed to the specified correlation functional.  If
  supplied, this takes precedence over additional arguments
  supplied via the \code{\dots} argument.}

  \item{useL1Median}{a logical indicating whether the
  \eqn{L_{1}}{L1} medians should be used as the centers of
  the data sets (defaults to \code{TRUE}).  If
  \code{FALSE}, the columnwise centers are used instead
  (columnwise means if \code{method} is \code{"pearson"}
  and columnwise medians otherwise).}

  \item{fast}{a logical indicating whether the faster
  modification of the algorithm should be used (defaults to
  \code{TRUE}).  See \dQuote{Details} for more
  information.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed to the
  specified correlation functional.}
}
\value{
  An object of class \code{"cca"} with the following
  components:

  \item{cor}{a numeric vector giving the canonical
  correlation measures.}

  \item{A}{a numeric matrix in which the columns contain
  the canonical vectors for \code{x}.}

  \item{B}{a numeric matrix in which the columns contain
  the canonical vectors for \code{y}.}

  \item{call}{the matched function call.}
}
\description{
  Perform canoncial correlation analysis via projection
  pursuit based on projections through the data points,
  with a focus on robust and nonparametric methods.
}
\details{
  First the candidate projection directions are defined for
  each data set from the respective center through each
  data point.  Then the basic algorithm scans all \eqn{n^2}
  possible combinations for the maximum correlation, where
  \eqn{n} is the number of observations.

  For larger values of \eqn{n}, the basic algorithm is
  computationally expensive, thus a faster modification is
  available as well.
}
\examples{
## generate data
library("mvtnorm")
set.seed(1234)  # for reproducibility
p <- 3
q <- 2
m <- p + q
sigma <- 0.5^t(sapply(1:m, function(i, j) abs(i-j), 1:m))
xy <- rmvnorm(100, sigma=sigma)
x <- xy[, 1:p]
y <- xy[, (p+1):m]

## Spearman correlation
ccaProj(x, y, method = "spearman")
ccaProj(x, y, method = "spearman", consistent = TRUE)

## Pearson correlation
ccaProj(x, y, method = "pearson")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{ccaGrid}} \code{\link{corFunctions}}
}
\keyword{multivariate}
\keyword{robust}

