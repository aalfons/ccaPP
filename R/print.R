# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' @S3method print cca
print.cca <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print maximum correlation
    cat("\nCanonical correlations:\n")
    print(x$cor, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print cvFolds
print.cvFolds <- function(x, ...) {
  # print general information
  cvText <- if(x$n == x$K) "Leave-one-out CV" else sprintf("%d-fold CV", x$K)
  if(x$R > 1) cvText <- paste(cvText, "with", x$R, "replications")
  cat(paste("\n", cvText, ":", sep=""))
  # print information on folds (add space between folds and subsets)
  subsets <- x$subsets
  if(x$R == 1) {
    cn <- if("grouping" %in% names(x)) "Group index" else "Index"
    nblanks <- 2
  } else {
    cn <- as.character(seq_len(x$R))
    nblanks <- 3
  }
  nblanks <- max(nchar(as.character(subsets[, 1]))-nchar(cn[1]), 0) + nblanks
  cn[1] <- paste(c(rep.int(" ", nblanks), cn[1]), collapse="")
  dimnames(subsets) <- list(Fold=x$which, cn)
  print(subsets, ...)
  # return object invisibly
  invisible(x)
}

#' @S3method print maxCor
print.maxCor <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print maximum correlation
    cat("\nMaximum correlation:\n")
    print(x$cor, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print permTest
print.permTest <- function(x, ...) {
    # print general statement
    cat("\nPermutation test for independence\n\n")
    # print maximum correlation and p-value
    cat(sprintf("r = %f, p-value = %f\n", x$cor0, x$pValue))
    # print number of random permuations
    cat(sprintf("R = %d random permuations\n", x$R))
    # print alternative hypothesis
    cat("Alternative hypothesis: true maximum correlation is not equal to 0\n")
    # return object invisibly
    invisible(x)
}

#' @S3method print sMaxCor
print.sMaxCor <- function(x, ...) {
  print.maxCor(x, ...)
  # print optimal tuning parameters
  if(!is.null(cv <- x$cv)) {
    haveGridX <- any(cv[, "lambdaX"] != 0)
    haveGridY <- any(cv[, "lambdaY"] != 0)
    if(haveGridX && haveGridY) cat("\nOptimal tuning parameters:\n")
    else cat("\nOptimal tuning parameter:\n")
    if(haveGridX) cat(sprintf("lambdaX: %f\n", x$lambdaX))
    if(haveGridY) cat(sprintf("lambdaY: %f\n", x$lambdaY))
  }
  # return object invisibly
  invisible(x)
}
