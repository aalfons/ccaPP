# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' (Robust) maximum correlation via alternating series of grid searches
#' 
#' Compute the maximum correlation between two data sets via projection pursuit 
#' based on alternating series of grid searches in two-dimensional subspaces of 
#' each data set, with a focus on robust and nonparametric methods.
#' 
#' \code{maxCorGrid} is based on alternating series of grid searches in 
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
#' Note that also the data sets are ordered according to the maximum average 
#' absolute correlation with the respective other data set to ensure symmetry 
#' of the algorithm.
#' 
#' \code{sMaxCorGrid} enforces sparsity by adding \eqn{L_{1}}{L1} penalties on 
#' the weighting vectors to the objective function (i.e., the correlation of 
#' the resulting projections).  The optimal combination of penalty paramters is 
#' thereby determined via cross-validation over a grid of values.  To increase 
#' the computational performance of this cross-validation procedure, parallel 
#' computing is implemented via package \pkg{parallel}.
#' 
#' @aliases print.maxCor
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param lambdaX,lambdaY  numeric vectors of non-negative values giving the 
#' penalty parameters for the weighting vectors of \code{x} and \code{y}, 
#' respectively.
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
#' @param standardize  a logical indicating whether the data should be 
#' (robustly) standardized.
#' @param fallback  logical indicating whether a fallback mode for robust 
#' standardization should be used.  If a correlation functional other than the 
#' Pearson correlation is maximized, the first attempt for standardizing the 
#' data is via median and MAD.  In the fallback mode, variables whose MADs are 
#' zero (e.g., dummy variables) are standardized via mean and standard 
#' deviation.  Note that if the Pearson correlation is maximized, 
#' standardization is always done via mean and standard deviation.
#' @param K,R,type,grouping   additional arguments for generating 
#' cross-validation folds (see \code{\link{cvFolds}}).
#' @param folds  an object of class \code{"cvFolds"} (as returned by 
#' \code{\link{cvFolds}}) defining the folds of the data for cross-validation.  
#' If supplied, this is preferred over the arguments for generating 
#' cross-validation folds.
#' @param nCores  a positive integer giving the number of processor cores to be 
#' used for parallel computing in cross-validation (the default is 1 for no 
#' parallelization).  If this is set to \code{NA}, all available processor 
#' cores are used.
#' @param cl  a \pkg{parallel} cluster for parallel computing to be used in 
#' cross-validation, as generated by \code{\link[parallel]{makeCluster}}.  If 
#' supplied, this is preferred over \code{nCores}.
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
#' @returnItem centerX  a numeric vector giving the center estimates used in 
#' standardization of \code{x}.
#' @returnItem centerY  a numeric vector giving the center estimates used in 
#' standardization of \code{y}.
#' @returnItem scaleX  a numeric vector giving the scale estimates used in 
#' standardization of \code{x}.
#' @returnItem scaleY  a numeric vector giving the scale estimates used in 
#' standardization of \code{y}.
#' @returnItem lambdaX  a numeric giving the optimal penalty parameter for the 
#' weighting vector of \code{x} (only \code{sMaxCorGrid}).
#' @returnItem lambdaY  a numeric giving the optimal penalty parameter for the 
#' weighting vector of \code{y} (only \code{sMaxCorGrid}).
#' @returnItem objective  a numeric giving the value of the objective function 
#' for the optimal combination of penalty parameters (only \code{sMaxCorGrid}).
#' @returnItem cv  a numeric matrix containing the (average) results from 
#' cross-validation for each combination of penalty parameters (only 
#' \code{sMaxCorGrid} and if a grid of values for the penalty parameters has 
#' been supplied).
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{maxCorProj}}, \code{\link{ccaGrid}}, 
#' \code{\link{corFunctions}}
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
                       control = list(...), nIterations = 10, nAlternate = 10, 
                       nGrid = 25, select = NULL, tol = 1e-06, 
                       standardize = TRUE, fallback = FALSE, seed = NULL, ...) {
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
                     algorithm="grid", ppControl=ppControl, 
                     standardize=standardize, fallback=fallback, 
                     seed=seed)
  maxCor$call <- matchedCall
  maxCor
}


#' @rdname maxCorGrid
#' @export

sMaxCorGrid <- function(x, y, lambdaX = 0, lambdaY = 0, 
                        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                        control = list(...), nIterations = 10, nAlternate = 10, 
                        nGrid = 25, select = NULL, tol = 1e-06, 
                        standardize = TRUE, fallback = FALSE, K = 5, R = 1, 
                        type = c("random", "consecutive", "interleaved"), 
                        grouping = NULL, folds = NULL, nCores = 1, cl = NULL, 
                        seed = NULL, ...) {
  ## initializations
  matchedCall <- match.call()
  ## define list of control arguments for algorithm
  lambdaX <- as.numeric(lambdaX)
  lambdaX[is.na(lambdaX)] <- formals()$lambdaX
  lambdaY <- as.numeric(lambdaY)
  lambdaY[is.na(lambdaY)] <- formals()$lambdaY
  nIterations <- as.integer(nIterations)
  nAlternate <- as.integer(nAlternate)
  nGrid <- as.integer(nGrid)
  tol <- as.numeric(tol)
  ppControl <- list(lambdaX=lambdaX, lambdaY=lambdaY, nIterations=nIterations, 
                    nAlternate=nAlternate, nGrid=nGrid, select=select, tol=tol)
  ## call workhorse function
  maxCor <- sMaxCorPP(x, y, method=method, corControl=control, 
                      algorithm="grid", ppControl=ppControl, 
                      standardize=standardize, fallback=fallback, K=K, 
                      R=R, type=type, grouping=grouping, folds=folds, 
                      nCores=nCores, cl=cl, seed=seed)
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
#' @param standardize  a logical indicating whether the data should be 
#' (robustly) standardized.
#' @param useL1Median  a logical indicating whether the \eqn{L_{1}}{L1} medians 
#' should be used as the centers of the data sets in standardization (defaults 
#' to \code{TRUE}).  If \code{FALSE}, the columnwise centers are used instead 
#' (columnwise means if \code{method} is \code{"pearson"} and columnwise 
#' medians otherwise).
#' @param fallback  logical indicating whether a fallback mode for robust 
#' standardization should be used.  If a correlation functional other than the 
#' Pearson correlation is maximized, the first attempt for standardizing the 
#' data is via median and MAD.  In the fallback mode, variables whose MADs are 
#' zero (e.g., dummy variables) are standardized via mean and standard 
#' deviation.  Note that if the Pearson correlation is maximized, 
#' standardization is always done via mean and standard deviation.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass maxCor
#' @returnItem cor  a numeric giving the maximum correlation estimate.
#' @returnItem a  numeric; the weighting vector for \code{x}.
#' @returnItem b  numeric; the weighting vector for \code{y}.
#' @returnItem centerX  a numeric vector giving the center estimates used in 
#' standardization of \code{x}.
#' @returnItem centerY  a numeric vector giving the center estimates used in 
#' standardization of \code{y}.
#' @returnItem scaleX  a numeric vector giving the scale estimates used in 
#' standardization of \code{x}.
#' @returnItem scaleY  a numeric vector giving the scale estimates used in 
#' standardization of \code{y}.
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
                       control = list(...), standardize = TRUE, 
                       useL1Median = TRUE, fallback = FALSE, ...) {
  ## initializations
  matchedCall <- match.call()
  ## define list of control arguments for algorithm
  ppControl <- list(useL1Median=isTRUE(useL1Median))
  ## call workhorse function
  maxCor <- maxCorPP(x, y, method=method, corControl=control, 
                     algorithm="proj", ppControl=ppControl, 
                     standardize=standardize, fallback=fallback)
  maxCor$call <- matchedCall
  maxCor
}


## workhorse function for maximum correlation
maxCorPP <- function(x, y, ...) {
  ## call workhorse function for canonical correlation analysis
  maxCor <- ccaPP(x, y, forceConsistency=FALSE, ...)
  ## modify object and return results
  maxCor <- list(cor=maxCor$cor, a=drop(maxCor$A), b=drop(maxCor$B), 
                 centerX=maxCor$centerX, centerY=maxCor$centerY, 
                 scaleX=maxCor$scaleX, scaleY=maxCor$scaleY)
  class(maxCor) <- "maxCor"
  maxCor
}


## workhorse function for sparse maximum correlation
#' @import parallel
sMaxCorPP <- function(x, y, 
                      method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                      corControl, algorithm = "grid", ppControl, 
                      standardize = TRUE, fallback = FALSE, folds = NULL, 
                      nCores = 1, cl = NULL, seed = NULL, ...) {
  ## initializations
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  if(nrow(y) != n) stop("'x' and 'y' must have the same number of observations")
  p <- ncol(x)
  q <- ncol(y)
  ## prepare the data and call C++ function
  if(n == 0 || p == 0 || q == 0) {
    # zero dimension
    maxCor <- list(cor=NA, a=numeric(), b=numeric(), lambdaX=NA, 
                   lambdaY=NA, objective=NA)
  } else {
    # get list of control arguments
    method <- match.arg(method)
    corControl <- getCorControl(method, corControl, forceConsistency=FALSE)
    standardize <- isTRUE(standardize)
    fallback <- isTRUE(fallback)
    if(!is.null(seed)) set.seed(seed)
    ppControl <- getPPControl(algorithm, ppControl, p=p, q=q)
    # check grid of values for the penalty parameters
    lambdaX <- ppControl$lambdaX
    lambdaY <- ppControl$lambdaY
    ppControl$lambdaX <- NULL
    ppControl$lambdaY <- NULL
    lambdaX <- if(p == 1) 0 else sort(unique(lambdaX))
    lambdaY <- if(q == 1) 0 else sort(unique(lambdaY))
    lambda <- as.matrix(expand.grid(lambdaX=lambdaX, lambdaY=lambdaY))
    ppControl$lambda <- lambda
    # check whether cross-validation over grid of tuning parameters is necessary
    useCV <- nrow(lambda) > 1
    if(useCV) {
      if(is.null(folds)) folds <- cvFolds(n, ...)
      # set up parallel computing if requested
      haveCl <- inherits(cl, "cluster")
      if(haveCl) haveNCores <- FALSE
      else {
        if(is.na(nCores)) nCores <- detectCores()  # use all available cores
        if(!is.numeric(nCores) || is.infinite(nCores) || nCores < 1) {
          nCores <- 1  # use default value
          warning("invalid value of 'nCores'; using default value")
        } else nCores <- as.integer(nCores)
        R <- folds$R
        nCores <- if(R == 1) min(nCores, folds$K) else min(nCores, R)
        haveNCores <- nCores > 1
      }
      # check whether parallel computing should be used
      useParallel <- haveNCores || haveCl
      # set up multicore or snow cluster if not supplied
      if(haveNCores) {
        if(.Platform$OS.type == "windows") {
          cl <- makePSOCKcluster(rep.int("localhost", nCores))
        } else cl <- makeForkCluster(nCores)
        on.exit(stopCluster(cl))
      }
      # perform cross-validation over grid of tuning parameters
      call <- call("sMaxCorPPFit", method=method, corControl=corControl, 
                   algorithm=algorithm, ppControl=ppControl, 
                   standardize=standardize, fallback=fallback)
      call[[1]] <- sMaxCorPPFit  # necessary to work with parallel computing
      corFun <- switch(method, spearman=corSpearman, kendall=corKendall, 
                       quadrant=corQuadrant, M=corM, pearson=corPearson)
      cv <- cvMaxCor(call, x, y, folds=folds, corFun=corFun, 
                     corArgs=corControl, cl=cl)
      # determine the optimal tuning parameters
      ppControl$lambda <- lambda[which.max(cv), , drop=FALSE]
    }
    # compute the final solution
    maxCor <- sMaxCorPPFit(x, y, method=method, corControl=corControl, 
                           algorithm=algorithm, ppControl=ppControl, 
                           standardize=standardize, fallback=fallback)
    maxCor$a <- drop(maxCor$a)
    maxCor$b <- drop(maxCor$b)
    if(useCV) maxCor$cv <- cbind(lambda, CV=cv)
  }
  ## assign class and return results
  class(maxCor) <- c("sMaxCor", "maxCor")
  maxCor
}

## simple wrapper for calling the C++ function
sMaxCorPPFit <- function(x, y, 
                         method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                         corControl, algorithm = "grid", ppControl, 
                         standardize = TRUE, fallback = FALSE) {
  # call C++ function
  maxCor <- .Call("R_sMaxCorPP", R_x=x, R_y=y, R_method=method, 
                  R_corControl=corControl, R_algorithm=algorithm, 
                  R_ppControl=ppControl, R_standardize=standardize, 
                  R_fallback=fallback, PACKAGE="ccaPP")
}
