# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
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

#' @export
print.permTest <- function(x, ...) {
    # print general statement
    cat("\nPermutation test for independence\n\n")
    # print maximum correlation and p-value
    cat(sprintf("p = %f, p-value = %f\n", x$cor, x$pval))
    # print alternative hypothesis
    cat("Alternative hypothesis: true maximum correlation is not equal to 0\n")
    # return object invisibly
    invisible(x)
}
