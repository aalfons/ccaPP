\name{ccaPP-package}
\alias{ccaPP-package}
\alias{ccaPP}
\docType{package}
\title{
(Robust) canonical correlation analysis via projection pursuit
}
\description{
Canonical correlation analysis and maximum correlation via 
    projection pursuit, as well as fast implementations of correlation
    estimators, with a focus on robust and nonparametric methods.
}
\details{
\tabular{ll}{
Package: \tab ccaPP\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2013-10-08\cr
Depends: \tab Rcpp (>= 0.9.10), parallel\cr
Imports: \tab RcppArmadillo (>= 0.2.36), pcaPP (>= 1.8-1), robustbase\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
Suggests: \tab mvtnorm\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
ccaGrid                 (Robust) CCA via alternating series of grid
                        searches
ccaMCD                  Robust CCA via alternating series of grid
                        searches based on the bivariate MCD
ccaPP-package           (Robust) canonical correlation analysis via
                        projection pursuit
ccaProj                 (Robust) CCA via projections through the data
                        points
corFunctions            Fast implementations of (robust) correlation
                        estimators
corMCD                  Robust correlation based on the bivariate
                        MCD-estimator
fastMAD                 Fast implementation of the median absolute
                        deviation
fastMedian              Fast implementation of the median
maxCorGrid              (Robust) maximum correlation via alternating
                        series of grid searches
maxCorMCD               Robust maximum correlation via alternating
                        series of grid searches based on the bivariate
                        MCD
maxCorProj              (Robust) maximum correlation via projections
                        through the data points
permTest                (Robust) permutation test for independence
}
}
\author{
Andreas Alfons [aut, cre],
  David Simcha [ctb]

Maintainer: Andreas Alfons <andreas.alfons@kuleuven.be>
}
\keyword{package}
