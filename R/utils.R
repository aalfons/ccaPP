# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

### transform canonical vectors back to original scale
#backtransform <- function(A, scale) {
#    apply(A, 2, 
#        function(a) {
#            sa <- a / scale       # divide by scale of corresponding variable
#            sa / sqrt(sum(sa^2))  # divide by norm
#        })
#}

## check if indices are within the limits
checkIndices <- function(indices, max) {
    indices <- as.integer(indices)
    indices[which(indices > 0 & indices <= max)]
}    

## call a function by either
# 1) simply evaluating a supplied function for the basic arguments if there are
#    no additional arguments in list format
# 2) evaluating a supplied function with 'do.call' if there are additional 
#    arguments in list format
doCall <- function(fun, ..., args = list()) {
  if(length(args) == 0) {
    fun(...)
  } else do.call(fun, c(list(...), args))
}

## call C++ to compute ranks of observations in a vector (for testing)
fastRank <- function(x) {
    x <- as.numeric(x)
    if(length(x) == 0) return(numeric())        # zero length vector
    .Call("R_rank", R_x=x, PACKAGE="ccaPP")   # call C++ function
}

## get list of control arguments for correlation function
getCorControl <- function(method, control, forceConsistency = TRUE) {
    if(method %in% c("spearman", "kendall", "quadrant")) {
        if(forceConsistency) out <- list(consistent=TRUE)
        else {
            # get default values (the three functions have the same arguments)
            out <- formals(corSpearman)[-(1:2)]
            # check supplied values
            if(is.list(control)) {
                if(!is.null(consistent <- control$consistent)) {
                    out$consistent <- isTRUE(consistent)
                }
            }
        }
    } else if(method == "M") {
        # get default values
        out <- formals(corM)[-(1:2)]
        choices <- eval(out$initial)
        out$initial <- choices[1]
        # check supplied values
        if(is.list(control)) {
            if(!is.null(prob <- control$prob)) {
                out$prob <- as.numeric(prob)
            }
            if(!is.null(initial <- control$initial)) {
                out$initial <- match.arg(initial, choices)
            }
            if(!is.null(tol <- control$tol)) {
                out$tol <- as.numeric(tol)
            }
        }
    } else out <- list()  # this case includes Pearson correlation
    # return list of control arguments
    out
}

## get list of control arguments for projection pursuit algorithm
getPPControl <- function(algorithm, control, p, q, seed = NULL) {
  # additional checks for grid search algorithm
  if(algorithm == "grid") {
    # check subset of variables to be used for determining the order of 
    # the variables from the respective other data set
    select <- control$select
    control$select <- NULL
    if(!is.null(select)) {
      if(is.list(select)) {
        # make sure select is a list with two index vectors and 
        # drop invalid indices from each vector
        select <- rep(select, length.out=2)
        select <- mapply(function(indices, max) {
          indices <- as.integer(indices)
          indices[which(indices > 0 & indices <= max)] - 1
        }, select, c(p, q))
        valid <- sapply(select, length) > 0
        # add the two index vectors to control object
        if(all(valid)) {
          control$selectX <- select[[1]]
          control$selectY <- select[[2]]
        } else select <- NULL
      } else {
        # check number of indices to sample
        select <- rep(as.integer(select), length.out=2)
        valid <- !is.na(select) & select > 0 & select < c(p, q)
        if(all(valid)) {
          # generate index vectors and add them to control object
          if(!is.null(seed)) set.seed(seed)
          control$selectX <- sample.int(p, select[1]) - 1
          control$selectY <- sample.int(q, select[2]) - 1
        } else select <- NULL
      }
    }
    if(is.null(select)) {
      control$selectX <- control$selectY <- integer()
    }
  }
  # return list of control arguments
  control
}

## L1 median (for testing)
l1Median <- function(x) {
    # initializations
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if(p == 0) return(numeric())       # no columns
    if(n == 0) return(rep.int(NA, p))  # no observations
    # call C++ function
    .Call("R_l1Median", R_x=x, PACKAGE="ccaPP")
}

### (robustly) standardize the data
#standardize <- function(x, robust = TRUE) {
#    if(robust) {
#        # with median and MAD
#        tmp <- apply(x, 2, function(v) unlist(fastMAD(v)))
#        center <- tmp[1,]  # column medians
#        x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
#        scale <- tmp[2,]  # column MADs
#        x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
#    } else {
#        # with mean and standard deviation
#        center <- colMeans(x)  # compute column means (faster than apply)
#        x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
#        f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
#        scale <- apply(x, 2, f)  # compute column scales with zero means
#        x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
#    }
#    # add attributes and return standardized data
#    attr(x, "center") <- center
#    attr(x, "scale") <- scale
#    x
#}
