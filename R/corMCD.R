# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @export
#' @import robustbase

corMCD <- function(x, y, alpha = 0.5) {
  covMcd(cbind(x, y), cor=TRUE, alpha=alpha)$cor[1, 2]
}

# # for testing
# corMcd <- function(x, y, alpha = 0.5) {
#   # initializations
#   x <- as.numeric(x)
#   y <- as.numeric(y)
#   alpha <- as.numeric(alpha)
#   # call C++ function
#   .Call("R_corMcd", R_x=x, R_y=y, R_alpha=alpha, PACKAGE="ccaPP")
# }
