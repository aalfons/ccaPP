# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## TODO: define class for results along with print method
#' @export
permTest <- function(x, y, R = 1000, algorithm = c("grid", "proj"), ..., 
        seed = NULL) {
    ## initializations
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if(nrow(y) != n) {
        stop("'x' and 'y' must have the same number of observations")
    }
    R <- rep(as.integer(R), length.out=1)
    if(is.na(R)) R <- formals()$R
    algorithm <- match.arg(algorithm)
    ccaFun <- switch(algorithm, grid=ccaGrid, proj=ccaProj)
    if(!is.null(seed)) set.seed(seed)  # set seed for reproducibility
    ## compute maximum correlation
    call <- as.call(list(ccaFun, ...))
    call$k <- 1  # compute only the first canonical correlation measure
    call$x <- x
    call$y <- y
    r <- eval(call)$cor
    ## permute rows of x and compute maximum correlation with y for each 
    ## permuted data set
    rPerm <- replicate(R, {
            call$x <- x[sample.int(n), , drop=FALSE]
            eval(call)$cor
        })
    ## compute p-value
    pval <- mean(r < rPerm)
    ## return results
    out <- list(pval=pval, cor=r)
    class(out) <- "permTest"
    out
}
