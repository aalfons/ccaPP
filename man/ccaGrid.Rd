\name{ccaGrid}
\alias{ccaGrid}
\alias{print.cca}
\title{(Robust) CCA via alternating series of grid searches}
\usage{
  ccaGrid(x, y, k = 1,
    method = c("spearman", "kendall", "quadrant", "M", "pearson"),
    control = list(...), nIterations = 10, nAlternate = 10,
    nGrid = 25, tol = 1e-06, ...)
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

  \item{nIterations}{an integer giving the maximum number
  of iterations.}

  \item{nAlternate}{an integer giving the maximum number of
  alternate grid searches in each iteration.}

  \item{nGrid}{an integer giving the number of equally
  spaced grid points on the unit circle to use in each grid
  search.}

  \item{tol}{a small positive numeric value to be used for
  determining convergence.}

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
  pursuit based on alternating series of grid searches in
  two-dimensional subspaces of each data set, with a focus
  on robust and nonparametric methods.
}
\details{
  The algorithm is based on alternating series of grid
  searches in two-dimensional subspaces of each data set.
  In each grid search, \code{nGrid} grid points on the unit
  circle in the corresponding plane are obtained, and the
  directions from the center to each of the grid points are
  examined.  In the first iteration, equispaced grid points
  in the interval \eqn{[-\pi/2, \pi/2)}{[-pi/2, pi/2)} are
  used.  In each subsequent iteration, the angles are
  halved such that the interval \eqn{[-\pi/4,
  \pi/4)}{[-pi/4, pi/4)} is used in the second iteration
  and so on.  If only one data set is multivariate, the
  algorithm simplifies to iterative grid searches in
  two-dimensional subspaces of the corresponding data set.
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
ccaGrid(x, y, method = "spearman")
ccaGrid(x, y, method = "spearman", consistent = TRUE)

## Pearson correlation
ccaGrid(x, y, method = "pearson")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{ccaProj}}, \code{\link{corFunctions}}
}
\keyword{multivariate}
\keyword{robust}

