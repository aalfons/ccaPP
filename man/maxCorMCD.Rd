\name{maxCorMCD}
\alias{maxCorMCD}
\title{Robust maximum correlation via alternating series of grid searches based on
the bivariate MCD}
\usage{
  maxCorMCD(x, y, alpha = 0.5, nIterations = 10,
    nAlternate = 10, nGrid = 25, tol = 1e-06,
    standardize = TRUE, fallback = FALSE)
}
\arguments{
  \item{x,y}{each can be a numeric vector, matrix or data
  frame.}

  \item{alpha}{a numeric value controlling the percentage
  of observations over which the determinant is minimized.
  Allowed values are between 0.5 and 1 and the default is
  0.5.}

  \item{nIterations}{an integer giving the maximum number
  of iterations.}

  \item{nAlternate}{an integer giving the maximum number of
  alternate series of grid searches in each iteration.}

  \item{nGrid}{an integer giving the number of equally
  spaced grid points on the unit circle to use in each grid
  search.}

  \item{tol}{a small positive numeric value to be used for
  determining convergence.}

  \item{fallback}{logical; the data are first robustly
  standardized via median and MAD.  This indicates whether
  standardization via mean and standard deviation should be
  performed as a fallback mode for variables whose MAD is
  zero (e.g., for dummy variables).}
}
\value{
  An object of class \code{"maxCor"} with the following
  components:

  \item{cor}{a numeric giving the maximum correlation
  estimate.}

  \item{a}{numeric; the weighting vector for \code{x}.}

  \item{b}{numeric; the weighting vector for \code{y}.}

  \item{call}{the matched function call.}
}
\description{
  Compute the maximum correlation between two data sets via
  projection pursuit based on alternating series of grid
  searches in two-dimensional subspaces of each data set.
  The correlations in those subspaces are computed with the
  bivariate MCD-estimator.
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

  The order of the variables in a series of grid searches
  for each of the data sets is determined by the average
  absolute correlations with the variables of the
  respective other data set.  This requires to compute the
  full \eqn{(p \times q)}{(p x q)} matrix of absolute
  correlations, where \eqn{p} denotes the number of
  variables of \code{x} and \eqn{q} the number of variables
  of \code{y}.
}
\note{
  Unlike \code{\link{maxCorGrid}}, \code{maxCorMCD} is a
  pure \R implementation of the alternate grid algorithm.
  For this reason, but mostly due to the computational cost
  of the MCD algorithm, \code{maxCorMCD} is much slower
  than \code{\link{maxCorGrid}}.
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

## maximum correlation
maxCorMCD(x, y, alpha = 0.75)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{maxCorGrid}}, \code{\link{ccaMCD}},
  \code{\link{corMCD}}
}
\keyword{multivariate}
\keyword{robust}

