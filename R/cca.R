# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @import Rcpp
#' @import RcppArmadillo
#' @import pcaPP
#' @useDynLib ccaPP

## CCA based on alternate grid searches
#' @export
ccaGrid <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), nIterations = 10, nAlternate = 10, nGrid = 25, 
        tol = 1e-06, ...) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    nIterations <- as.integer(nIterations)
    nAlternate <- as.integer(nAlternate)
    nGrid <- as.integer(nGrid)
    tol <- as.numeric(tol)
    ppControl <- list(nIterations=nIterations, nAlternate=nAlternate, 
        nGrid=nGrid, tol=tol)
    ## call workhorse function
    cca <- ccaPP(x, y, k, method=method, corControl=control, 
        algorithm="grid", ppControl=ppControl)
    cca$call <- matchedCall
    cca
}

## CCA based on projections through the data points
#' @export
ccaProj <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), useL1Median = TRUE, fast = FALSE, ..., 
        seed = NULL) {
    ## initializations
    matchedCall <- match.call()
    ## define list of control arguments for algorithm
    ppControl <- list(useL1Median=isTRUE(useL1Median), fast=isTRUE(fast))
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
        if(k > 1) {
            k <- 1
            warning("currently only implemented for k = 1")
        }
        # check method and get list of control arguments
        method <- match.arg(method)
        corControl <- getCorControl(method, corControl)
        # get initial random starting indices for fast algorithm based on 
        # projections through data points
        if(algorithm == "proj" && ppControl$fast && p > 1 && q > 1) {
            if(!is.null(seed)) set.seed(seed)
            ppControl$initial <- sample.int(n, 2) - 1
        }
        # standardize the data
        if(method == "pearson") {
            x <- standardize(x, robust=FALSE)
            y <- standardize(y, robust=FALSE)
        } else {
            x <- standardize(x, robust=TRUE)
            y <- standardize(y, robust=TRUE)
        }
        # call C++ function
        cca <- .Call("R_ccaPP", R_x=x, R_y=y, R_k=k, R_method=method, 
            R_corControl=corControl, R_algorithm=algorithm, 
            R_ppControl=ppControl, PACKAGE="ccaPP")
        # transform canonical vectors back to original scale
        cca$A <- backtransform(cca$A, attr(x, "scale"))
        cca$B <- backtransform(cca$B, attr(y, "scale"))
}
    class(cca) <- "cca"
    cca
}
