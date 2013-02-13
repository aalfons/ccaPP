# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## perform cross-validation
cvMaxCor <- function(call, x, y, folds, corFun = corSpearman, corArgs = list(), 
                     envir = parent.frame(), cl = NULL) {
  # initializations
  R <- folds$R
  useParallel <- !is.null(cl)
  # obtain list of predictions for all replications
  if(R == 1) {
    s <- getIndices(folds)
    if(useParallel) xy <- pcvXY(s, call=call, x=x, y=y, envir=envir, cl=cl)
    else xy <- cvXY(s, call=call, x=x, y=y, envir=envir)
    corXY(xy, corFun=corFun, corArgs=corArgs)
  } else {
    cvFun <- function(r) {
      s <- getIndices(folds, r)
      xy <- cvXY(s, call=call, x=x, y=y, envir=envir)
      corXY(xy, corFun=corFun, corArgs=corArgs)
    }
    if(useParallel) {
      cv <- parLapply(cl, seq_len(R), cvFun)
      cv <- do.call(cbind, cv)
    }
    else cv <- sapply(seq_len(R), cvFun)
    apply(cv, 1, mean)
  }
}

# compute the weighting vectors for the training data and obtain linear 
# combinations of the test data
rsXY <- function(i, call, x, y, envir) {
  # plug training data into function call
  call$x <- x[-i, , drop=FALSE]
  call$y <- y[-i, , drop=FALSE]
  # evaluate function call in supplied environment to make sure 
  # that other arguments are evaluated correctly
  fit <- eval(call, envir)
  # obtain linear combinations for test data
  cbind(x=x[i, , drop=FALSE] %*% fit$a, y=y[i, , drop=FALSE] %*% fit$b)
}

# compute linear combinations for each left-out block in one replication of 
# cross validation
cvXY <- function(folds, call, x, y, envir) {
  tmp <- lapply(folds, rsXY, call=call, x=x, y=y, envir=envir)
  # instead of collecting the results from the folds in the original order 
  # of the observations, they are simply stacked on top of each other
  do.call(rbind, tmp)
}

# compute linear combinations for each left-out block in one replication of 
# cross validation via parallel computing
pcvXY <- function(folds, call, x, y, envir, cl) {
  # parLapply() already has an argument 'x', so don't use the argument names 
  # for the data or it will throw an error
  tmp <- parLapply(cl, folds, rsXY, call, x, y, envir=envir)
  # instead of collecting the results from the folds in the original order 
  # of the observations, they are simply stacked on top of each other
  do.call(rbind, tmp)
}

# compute the correlation between the linear combinations obtained via 
# cross-validation
corXY <- function(xy, corFun = corSpearman, corArgs = list()) {
  # number of linear combinations of each data set
  p <- ncol(xy) / 2
  # split the data into the linear combinations of x and y, respectively, and 
  # compute the correlation between each pair of linear combinations
  if(p == 1) doCall(corFun, xy[, 1], xy[, 2], args=corArgs)
  else {
    sapply(seq_len(p), function(j) {
      doCall(corFun, xy[, j], xy[, p+j], args=corArgs)
    })
  }
}


#' Cross-validation folds
#' 
#' Split observations or groups of observations into \eqn{K} folds to be used 
#' for (repeated) \eqn{K}-fold cross-validation.  \eqn{K} should thereby be 
#' chosen such that all folds are of approximately equal size.
#' 
#' @aliases print.cvFolds
#' 
#' @param n  an integer giving the number of observations to be split into 
#' folds.  This is ignored if \code{grouping} is supplied in order to split 
#' groups of observations into folds.
#' @param K  an integer giving the number of folds into which the observations 
#' should be split (the default is five).  Setting \code{K} equal to the number 
#' of observations or groups yields leave-one-out cross-validation.
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.  This is ignored for for leave-one-out 
#' cross-validation and other non-random splits of the data.
#' @param type  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' @param grouping  a factor specifying groups of observations.  If supplied, 
#' the data are split according to the groups rather than individual 
#' observations such that all observations within a group belong to the same 
#' fold.
#' 
#' @returnClass cvFolds
#' @returnItem n  an integer giving the number of observations or groups.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem subsets  an integer matrix in which each column contains a 
#' permutation of the indices of the observations or groups.
#' @returnItem which  an integer vector giving the fold for each permuted 
#' observation or group.
#' @returnItem grouping  a list giving the indices of the observations 
#' belonging to each group.  This is only returned if a grouping factor 
#' has been supplied.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{sMaxCorGrid}}
#' 
#' @examples 
#' set.seed(1234)  # set seed for reproducibility
#' cvFolds(20, K = 5)
#' cvFolds(20, K = 5, R = 10)
#' 
#' @keywords utilities
#' 
#' @export 

cvFolds <- function(n, K = 5, R = 1, 
                    type = c("random", "consecutive", "interleaved"), 
                    grouping = NULL) {
  # check arguments
  n <- if(is.null(grouping)) round(rep(n, length.out=1)) else nlevels(grouping)
  if(!isTRUE(n > 0)) stop("'n' must be positive")
  if(!isTRUE(K <= n)) stop(sprintf("'K' must be smaller or equal to %d", n))
  if(K == n) type <- "leave-one-out"
  else type <- match.arg(type)
  # obtain CV folds
  if(type == "random") {
    # random K-fold splits with R replications
    subsets <- replicate(R, sample.int(n))
  } else {
    # leave-one-out CV or non-random splits, replication not meaningful
    R <- 1
    subsets <- as.matrix(seq_len(n))
  }
  which <- as.factor(rep(seq_len(K), length.out=n))
  if(type == "consecutive") which <- rep.int(seq_len(K), tabulate(which))
  # construct and return object
  folds <- list(n=n, K=K, R=R, subsets=subsets, which=which)
  if(!is.null(grouping)) folds$grouping <- split(seq_along(grouping), grouping)
  class(folds) <- "cvFolds"
  folds
}

# retrieve list of indices for r-th replication
getIndices <- function(x, r = 1) {
  # split permuted items according to the folds
  subsets <- split(x$subsets[, r], x$which)
  # in case of grouped data, the list contains the group indices in each CV 
  # fold, so the indices of the respective observations need to be extracted
  if(!is.null(grouping <- x$grouping)) 
    subsets <- lapply(subsets, function(s) unlist(grouping[s], use.names=FALSE))
  # return list of indices for CV folds
  names(subsets) <- NULL
  subsets
}
