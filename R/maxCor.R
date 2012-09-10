# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' (Robust) maximum correlation via alternating series of grid searches
#' 
#' Compute the maximum correlation between two data sets via projection pursuit 
#' based on alternating series of grid searches in two-dimensional subspaces of 
#' each data set, with a focus on robust and nonparametric methods.
#' 
#' The algorithm is based on alternating series of grid searches in 
#' two-dimensional subspaces of each data set.  In each grid search, 
#' \code{nGrid} grid points on the unit circle in the corresponding plane are 
#' obtained, and the directions from the center to each of the grid points are 
#' examined.  In the first iteration, equispaced grid points in the interval 
#' \eqn{[-\pi/2, \pi/2)}{[-pi/2, pi/2)} are used.  In each subsequent 
#' iteration, the angles are halved such that the interval 
#' \eqn{[-\pi/4, \pi/4)}{[-pi/4, pi/4)} is used in the second iteration and so 
#' on.  If only one data set is multivariate, the algorithm simplifies 
#' to iterative grid searches in two-dimensional subspaces of the corresponding 
#' data set.
#' 
#' In the basic algorithm, the order of the variables in a series of grid 
#' searches for each of the data sets is determined by the average absolute 
#' correlations with the variables of the respective other data set.  Since 
#' this requires to compute the full \eqn{(p \times q)}{(p x q)} matrix of 
#' absolute correlations, where \eqn{p} denotes the number of variables of 
#' \code{x} and \eqn{q} the number of variables of \code{y}, a faster 
#' modification is available as well.  In this modification, the average 
#' absolute correlations are computed over only a subset of the variables of 
#' the respective other data set.  It is thereby possible to use randomly 
#' selected subsets of variables, or to specify the subsets of variables 
#' directly.
#' 
#' @aliases print.maxCor
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param method  a character string specifying the correlation functional to 
#' maximize.  Possible values are \code{"spearman"} for the Spearman 
#' correlation, \code{"kendall"} for the Kendall correlation, \code{"quadrant"} 
#' for the quadrant correlation, \code{"M"} for the correlation based on a 
#' bivariate M-estimator of location and scatter with a Huber loss function, or 
#' \code{"pearson"} for the classical Pearson correlation (see 
#' \code{\link{corFunctions}}).
#' @param control  a list of additional arguments to be passed to the specified 
#' correlation functional.  If supplied, this takes precedence over additional 
#' arguments supplied via the \code{\dots} argument.
#' @param nIterations  an integer giving the maximum number of iterations.
#' @param nAlternate  an integer giving the maximum number of alternate series 
#' of grid searches in each iteration.
#' @param nGrid  an integer giving the number of equally spaced grid points on 
#' the unit circle to use in each grid search.
#' @param select  optional; either an integer vector of length two or a list 
#' containing two index vectors.  In the first case, the first integer gives 
#' the number of variables of \code{x} to be randomly selected for determining 
#' the order of the variables of \code{y} in the corresponding series of grid 
#' searches, and vice versa for the second integer.  In the latter case, the 
#' first list element gives the indices of the variables of \code{x} to be used 
#' for determining the order of the variables of \code{y}, and vice versa for 
#' the second integer (see \dQuote{Details}).
#' @param tol  a small positive numeric value to be used for determining 
#' convergence.
#' @param fallback  logical; if a correlation functional other than the 
#' Pearson correlation is maximized, the data are first robustly standardized 
#' via median and MAD.  This indicates whether standardization via mean and 
#' standard deviation should be performed as a fallback mode for variables 
#' whose MAD is zero (e.g., for dummy variables).  Note that if the Pearson 
#' correlation is maximized, the data are always standardized via mean and 
#' standard deviation.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).  This is only used if \code{select} specifies 
#' the numbers of variables of each data set to be randomly selected for 
#' determining the order of the variables of the respective other data set.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass maxCor
#' @returnItem cor  a numeric giving the maximum correlation estimate.
#' @returnItem a  numeric; the weighting vector for \code{x}.
#' @returnItem b  numeric; the weighting vector for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{maxCorProj}}, \code{\link{ccaGrid}}, 
#' \code{\link{corFunctions}}, 
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' p <- 3
#' q <- 2
#' m <- p + q
#' sigma <- 0.5^t(sapply(1:m, function(i, j) abs(i-j), 1:m))
#' xy <- rmvnorm(100, sigma=sigma)
#' x <- xy[, 1:p]
#' y <- xy[, (p+1):m]
#' 
#' ## Spearman correlation
#' maxCorGrid(x, y, method = "spearman")
#' maxCorGrid(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' maxCorGrid(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib ccaPP
#' @export

maxCorGrid <- function(x, y, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), nIterations = 10, nAlternate = 10, nGrid = 25, 
        select = NULL, tol = 1e-06, fallback = FALSE, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    nIterations <- as.integer(nIterations)
    nAlternate <- as.integer(nAlternate)
    nGrid <- as.integer(nGrid)
    tol <- as.numeric(tol)
    ppControl <- list(nIterations=nIterations, nAlternate=nAlternate, 
        nGrid=nGrid, select=select, tol=tol)
    ## call workhorse function
    maxCor <- maxCorPP(x, y, method=method, corControl=control, 
        algorithm="grid", ppControl=ppControl, fallback=fallback, 
        seed=seed)
    maxCor$call <- matchedCall
    maxCor
}


#' (Robust) maximum correlation via projections through the data points
#' 
#' Compute the maximum correlation between two data sets via projection pursuit 
#' based on projections through the data points, with a focus on robust and 
#' nonparametric methods.
#' 
#' First the candidate projection directions are defined for each data set 
#' from the respective center through each data point.  Then the algorithm 
#' scans all \eqn{n^2} possible combinations for the maximum correlation, 
#' where \eqn{n} is the number of observations.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param method  a character string specifying the correlation functional to 
#' maximize.  Possible values are \code{"spearman"} for the Spearman 
#' correlation, \code{"kendall"} for the Kendall correlation, \code{"quadrant"} 
#' for the quadrant correlation, \code{"M"} for the correlation based on a 
#' bivariate M-estimator of location and scatter with a Huber loss function, or 
#' \code{"pearson"} for the classical Pearson correlation (see 
#' \code{\link{corFunctions}}).
#' @param control  a list of additional arguments to be passed to the specified 
#' correlation functional.  If supplied, this takes precedence over additional 
#' arguments supplied via the \code{\dots} argument.
#' @param useL1Median  a logical indicating whether the \eqn{L_{1}}{L1} medians 
#' should be used as the centers of the data sets (defaults to 
#' \code{TRUE}).  If \code{FALSE}, the columnwise centers are used instead 
#' (columnwise means if \code{method} is \code{"pearson"} and columnwise 
#' medians otherwise).
#' @param fallback  logical; if a correlation functional other than the 
#' Pearson correlation is maximized, the data are first robustly standardized 
#' via median and MAD.  This indicates whether standardization via mean and 
#' standard deviation should be performed as a fallback mode for variables 
#' whose MAD is zero (e.g., for dummy variables).  Note that if the Pearson 
#' correlation is maximized, the data are always standardized via mean and 
#' standard deviation.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass maxCor
#' @returnItem cor  a numeric giving the maximum correlation estimate.
#' @returnItem a  numeric; the weighting vector for \code{x}.
#' @returnItem b  numeric; the weighting vector for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{maxCorGrid}}, \code{\link{ccaProj}}, 
#' \code{\link{corFunctions}}, 
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' p <- 3
#' q <- 2
#' m <- p + q
#' sigma <- 0.5^t(sapply(1:m, function(i, j) abs(i-j), 1:m))
#' xy <- rmvnorm(100, sigma=sigma)
#' x <- xy[, 1:p]
#' y <- xy[, (p+1):m]
#' 
#' ## Spearman correlation
#' maxCorProj(x, y, method = "spearman")
#' maxCorProj(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' maxCorProj(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @import pcaPP
#' @useDynLib ccaPP
#' @export

maxCorProj <- function(x, y, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), useL1Median = TRUE, fallback = FALSE, ...) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    ppControl <- list(useL1Median=isTRUE(useL1Median))
    ## call workhorse function
    maxCor <- maxCorPP(x, y, method=method, corControl=control, 
        algorithm="proj", ppControl=ppControl, fallback=fallback)
    maxCor$call <- matchedCall
    maxCor
}


## workhorse function
maxCorPP <- function(x, y, ...) {
    ## call workhorse function for canonical correlation analysis
    maxCor <- ccaPP(x, y, forceConsistency=FALSE, ...)
    ## modify object and return results
    maxCor <- list(cor=maxCor$cor, a=drop(maxCor$A), b=drop(maxCor$B))
    class(maxCor) <- "maxCor"
    maxCor
}