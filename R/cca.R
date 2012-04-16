# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' (Robust) CCA via alternating series of grid searches
#' 
#' Perform canoncial correlation analysis via projection pursuit based on 
#' alternating series of grid searches in two-dimensional subspaces of each 
#' data set, with a focus on robust and nonparametric methods.
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
#' @aliases print.cca
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
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
#' @param nAlternate  an integer giving the maximum number of alternate grid 
#' searches in each iteration.
#' @param nGrid  an integer giving the number of equally spaced grid points on 
#' the unit circle to use in each grid search.
#' @param tol  a small positive numeric value to be used for determining 
#' convergence.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass cca
#' @returnItem cor  a numeric vector giving the canonical correlation 
#' measures.
#' @returnItem A  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{x}.
#' @returnItem B  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaProj}}, \code{\link{corFunctions}}
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
#' ccaGrid(x, y, method = "spearman")
#' ccaGrid(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' ccaGrid(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @import pcaPP
#' @useDynLib ccaPP
#' @export

ccaGrid <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), nIterations = 10, nAlternate = 10, nGrid = 25, 
        initial = NULL, tol = 1e-06, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    nIterations <- as.integer(nIterations)
    nAlternate <- as.integer(nAlternate)
    nGrid <- as.integer(nGrid)
    tol <- as.numeric(tol)
    ppControl <- list(nIterations=nIterations, nAlternate=nAlternate, 
        nGrid=nGrid, initial=initial, tol=tol)
    ## call workhorse function
    cca <- ccaPP(x, y, k, method=method, corControl=control, 
        algorithm="grid", ppControl=ppControl)
    cca$call <- matchedCall
    cca
}


#' (Robust) CCA via projections through the data points
#' 
#' Perform canoncial correlation analysis via projection pursuit based on 
#' projections through the data points, with a focus on robust and 
#' nonparametric methods.
#' 
#' First the candidate projection directions are defined for each data set 
#' from the respective center through each data point.  Then the basic 
#' algorithm scans all \eqn{n^2} possible combinations for the maximum 
#' correlation, where \eqn{n} is the number of observations.
#' 
#' For larger values of \eqn{n}, the basic algorithm is computationally 
#' expensive, thus a faster modification is available as well.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
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
#' @param fast  a logical indicating whether the faster modification of the 
#' algorithm should be used (defaults to \code{TRUE}).  See \dQuote{Details} 
#' for more information.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.
#' 
#' @returnClass cca
#' @returnItem cor  a numeric vector giving the canonical correlation 
#' measures.
#' @returnItem A  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{x}.
#' @returnItem B  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaGrid}} \code{\link{corFunctions}}
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
#' ccaProj(x, y, method = "spearman")
#' ccaProj(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' ccaProj(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @import Rcpp
#' @import RcppArmadillo
#' @import pcaPP
#' @useDynLib ccaPP
#' @export

ccaProj <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), useL1Median = TRUE, fast = FALSE, 
        nIterations = 10, initial = NULL, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    nIterations <- as.integer(nIterations)
    ## define list of control arguments for algorithm
    ppControl <- list(useL1Median=isTRUE(useL1Median), fast=isTRUE(fast), 
        nIterations=nIterations, initial=initial)
    ## call workhorse function
    cca <- ccaPP(x, y, k, method=method, corControl=control, algorithm="proj", 
        ppControl=ppControl, seed=seed)
    cca$call <- matchedCall
    cca
}


## workhorse function
ccaPP <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        corControl, algorithm = c("grid", "proj"), ppControl, 
        seed = NULL) {
    ## initializations
    matchedCall <- match.call()
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if(nrow(y) != n) {
        stop("'x' and 'y' must have the same number of observations")
    }
    p <- ncol(x)
    q <- ncol(y)
    # check number of canonical variables to compute
    k <- rep(as.integer(k), length.out=1)
    if(is.na(k) || k < 0) k <- formals()$k
    k <- min(k, p, q)
    ## prepare the data and call C++ function
    if(n == 0 || p == 0 || q == 0 || k == 0) {
        # zero dimension
        A <- B <- matrix(numeric(), 0, 0)
        cca <- list(cor=NA, A=A, B=B)
    } else {
        # check method and get list of control arguments
        method <- match.arg(method)
        corControl <- getCorControl(method, corControl)
        if(algorithm == "grid") {
            # check subset of variables to be used for determining the order of 
            # the variables from the respective other data set
            initial <- ppControl$initial
            ppControl$initial <- NULL
            if(!is.null(initial)) {
                if(is.list(initial)) {
                    # make sure initial is a list with two index vectors and 
                    # drop invalid indices from each vector
                    initial <- rep(initial, length.out=2)
                    initial <- mapply(function(indices, max) {
                            indices <- as.integer(indices)
                            indices[which(indices > 0 & indices <= max)] - 1
                        }, initial, c(p, q))
                    valid <- sapply(initial, length) > 0
                    # add the two index vectors to control object
                    if(all(valid)) {
                        ppControl$initialX <- initial[[1]]
                        ppControl$initialY <- initial[[2]]
                    } else initial <- NULL
                } else {
                    # check number of indices to sample
                    initial <- rep(as.integer(initial), length.out=2)
                    valid <- !is.na(initial) & initial > 0 & initial < c(p, q)
                    if(all(valid)) {
                        # generate index vectors and add them to control object
                        if(!is.null(seed)) set.seed(seed)
                        ppControl$initialX <- sample.int(p, initial[1]) - 1
                        ppControl$initialY <- sample.int(q, initial[2]) - 1
                    } else initial <- NULL
                }
            }
            if(is.null(initial)) {
                ppControl$initialX <- ppControl$initialY <- integer()
            }
        } else if(algorithm == "proj") {
            # get initial random starting indices for fast algorithm based on 
            # projections through data points
            if(ppControl$fast && p > 1 && q > 1) {
                if(is.null(ppControl$initial)) {
                    if(!is.null(seed)) set.seed(seed)
                    ppControl$initial <- sample.int(n, 2) - 1
                } else {
                    initial <- rep(as.integer(ppControl$initial), length.out=2)
                    # if one of the starting points has invalid value, replace 
                    # it with random starting point
                    invalid <- is.na(initial) | initial < 1 | initial > n
                    if(any(invalid)) {
                        if(!is.null(seed)) set.seed(seed)
                        if(invalid[1]) initial[1] <- sample.int(n, 1)
                        if(invalid[2]) initial[2] <- sample.int(n, 1)
                    }
                    ppControl$initial <- initial - 1
                }
            } else ppControl$initial <- integer()
        }
#        # standardize the data
#        if(method == "pearson") {
#            x <- standardize(x, robust=FALSE)
#            y <- standardize(y, robust=FALSE)
#        } else {
#            x <- standardize(x, robust=TRUE)
#            y <- standardize(y, robust=TRUE)
#        }
        # call C++ function
        cca <- .Call("R_ccaPP", R_x=x, R_y=y, R_k=k, R_method=method, 
            R_corControl=corControl, R_algorithm=algorithm, 
            R_ppControl=ppControl, PACKAGE="ccaPP")
#        # transform canonical vectors back to original scale
#        cca$A <- backtransform(cca$A, attr(x, "scale"))
#        cca$B <- backtransform(cca$B, attr(y, "scale"))
}
    class(cca) <- "cca"
    cca
}
