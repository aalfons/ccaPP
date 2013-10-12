\name{corMCD}
\alias{corMCD}
\title{Robust correlation based on the bivariate MCD-estimator}
\usage{
  corMCD(x, y, alpha = 0.5, seed = NULL)
}
\arguments{
  \item{x,y}{numeric vectors.}

  \item{alpha}{a numeric value controlling the percentage
  of observations over which the determinant is minimized.
  Allowed values are between 0.5 and 1 and the default is
  0.5.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}
}
\value{
  The correlation estimate.
}
\description{
  Estimate the correlation of two vectors based on the
  bivariate MCD-estimator of location and scatter.
}
\note{
  This is a simple wrapper function calling
  \code{\link[robustbase]{covMcd}} from package
  \pkg{robustbase}.
}
\examples{
## generate data
library("mvtnorm")
set.seed(1234)  # for reproducibility
sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
xy <- rmvnorm(100, sigma=sigma)
x <- xy[, 1]
y <- xy[, 2]

## compute correlation
corMCD(x, y, alpha = 0.75)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{ccaMCD}}, \code{\link[robustbase]{covMcd}}
}
\keyword{multivariate}
\keyword{robust}

