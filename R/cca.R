# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## CCA based on alternate grid searches
ccaGrid <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), nIterations = 10, nAlternate = 10, nGrid = 25, 
        tol = 1e-06, ...) {
    ## define list of control arguments for algorithm
    nIterations <- as.integer(nIterations)
    nAlternate <- as.integer(nAlternate)
    nGrid <- as.integer(nGrid)
    tol <- as.numeric(tol)
    ppControl <- list(nIterations=nIterations, nAlternate=nAlternate, 
        nGrid=nGrid, tol=tol)
    ## call workhorse function
    ccaPP(x, y, k, method=method, corControl=control, 
        algorithm="grid", ppControl=ppControl)
}

## CCA based on alternate grid searches
ccaProj <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        control = list(...), useL1Median = TRUE, fast = FALSE, ...) {
    ## define list of control arguments for algorithm
    algorithm <- if(isTRUE(fast)) "fast" else "naive"
    ppControl <- list(useL1Median=isTRUE(useL1Median))
    ## call workhorse function
    ccaPP(x, y, k, method=method, corControl=control, 
        algorithm=algorithm, ppControl=ppControl)
}

## workhorse function
ccaPP <- function(x, y, k = 1, 
        method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
        corControl, algorithm = c("grid", "naive", "fast"), ppControl) {
    ## initializations
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if(nrow(y) != n) {
        stop("'x' and 'y' must have the same number of observations")
    }
    p <- ncol(x)
    q <- ncol(y)
    if(n == 0 || p == 0 || q == 0) {
        # zero dimension
        A <- B <- matrix(numeric(), 0, 0)
        cca <- list(cor=NA, A=A, B=B)
    } else {
        # check number of canonical variables to compute
        k <- rep(as.integer(k), length.out=1)
        if(is.na(k)) k <- formals()$k
        k <- min(k, p, q)
        # check method and get list of control arguments
        method <- match.arg(method)
        corControl <- getCorControl(method, corControl)
        # standardize the data
        if(method == "pearson") {
            x <- standardize(x, robust=FALSE)
            y <- standardize(y, robust=FALSE)
        } else {
            x <- standardize(x, robust=TRUE)
            y <- standardize(y, robust=TRUE)
        }
        stop("not implemented yet")
    }
}
