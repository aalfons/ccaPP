# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Robust correlation based on the bivariate MCD-estimator
#' 
#' Estimate the correlation of two vectors based on the bivariate MCD-estimator 
#' of location and scatter.
#' 
#' @param x,y  numeric vectors.
#' @param alpha  a numeric value controlling the percentage of observations 
#' over which the determinant is minimized.  Allowed values are between 0.5 
#' and 1 and the default is 0.5.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' 
#' @return The correlation estimate.
#' 
#' @note This is a simple wrapper function calling 
#' \code{\link[robustbase]{covMcd}} from package \pkg{robustbase}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaMCD}}, \code{\link[robustbase]{covMcd}}
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
#' xy <- rmvnorm(100, sigma=sigma)
#' x <- xy[, 1]
#' y <- xy[, 2]
#' 
#' ## compute correlation
#' corMCD(x, y, alpha = 0.75)
#' 
#' @keywords multivariate robust
#' 
#' @export
#' @import robustbase

corMCD <- function(x, y, alpha = 0.5, seed = NULL) {
  xy <- cbind(as.numeric(x), as.numeric(y))
  covMcd(xy, cor=TRUE, alpha=alpha, seed=seed)$cor[1, 2]
}


#' Robust maximum correlation via alternating series of grid searches based on 
#' the bivariate MCD
#' 
#' Compute the maximum correlation between two data sets via projection pursuit 
#' based on alternating series of grid searches in two-dimensional subspaces of 
#' each data set.  The correlations in those subspaces are computed with the 
#' bivariate MCD-estimator.
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
#' The order of the variables in a series of grid searches for each of the data 
#' sets is determined by the average absolute correlations with the variables 
#' of the respective other data set.  This requires to compute the full 
#' \eqn{(p \times q)}{(p x q)} matrix of absolute correlations, where \eqn{p} 
#' denotes the number of variables of \code{x} and \eqn{q} the number of 
#' variables of \code{y}.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param alpha  a numeric value controlling the percentage of observations 
#' over which the determinant is minimized.  Allowed values are between 0.5 
#' and 1 and the default is 0.5.
#' @param nIterations  an integer giving the maximum number of iterations.
#' @param nAlternate  an integer giving the maximum number of alternate series 
#' of grid searches in each iteration.
#' @param nGrid  an integer giving the number of equally spaced grid points on 
#' the unit circle to use in each grid search.
#' @param tol  a small positive numeric value to be used for determining 
#' convergence.
#' @param fallback  logical; the data are first robustly standardized via 
#' median and MAD.  This indicates whether standardization via mean and 
#' standard deviation should be performed as a fallback mode for variables 
#' whose MAD is zero (e.g., for dummy variables).
#' 
#' @returnClass maxCor
#' @returnItem cor  a numeric giving the maximum correlation estimate.
#' @returnItem a  numeric; the weighting vector for \code{x}.
#' @returnItem b  numeric; the weighting vector for \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @note Unlike \code{\link{maxCorGrid}}, \code{maxCorMCD} is a pure \R 
#' implementation of the alternate grid algorithm.  For this reason, but mostly 
#' due to the computational cost of the MCD algorithm, \code{maxCorMCD} is much 
#' slower than \code{\link{maxCorGrid}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{maxCorGrid}}, \code{\link{ccaMCD}}, \code{\link{corMCD}}
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
#' ## maximum correlation
#' maxCorMCD(x, y, alpha = 0.75)
#' 
#' @keywords multivariate robust
#' 
#' @export

maxCorMCD <- function(x, y, alpha = 0.5, nIterations = 10, 
                      nAlternate = 10, nGrid = 25, tol = 1e-06, 
                      standardize = TRUE, fallback = FALSE) {
  ## initializations
  matchedCall <- match.call()
  x <- as.matrix(x)
  y <- as.matrix(y)
  standardize <- isTRUE(standardize)
  fallback <- isTRUE(fallback)
  ## check control arguments for PP algorithm
  nIterations <- as.integer(nIterations)
  nAlternate <- as.integer(nAlternate)
  nGrid <- as.integer(nGrid)
  tol <- as.numeric(tol)
  ## check control arguments for correlation based on bivariate MCD-estimator
  alpha <- as.numeric(alpha)
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) runif(1)
  seed <- .Random.seed
  on.exit(set.seed(seed))  # reset seed of the random number generator
  ## robustly standardize data if requested
  if(standardize) {
    xs <- robStandardize(x, fallback=fallback)
    ys <- robStandardize(y, fallback=fallback)
    scaleX <- attr(xs, "scale")
    scaleY <- attr(ys, "scale")
    ## return results
    maxCor <- maxCorPPMCD(xs, ys, alpha=alpha, nIterations=nIterations, 
                          nAlternate=nAlternate, nGrid=nGrid, tol=tol, 
                          seed=seed)
    maxCor <- list(cor=maxCor$cor, a=backtransform(maxCor$a, scaleX), 
                   b=backtransform(maxCor$b, scaleY))
  } else {
    maxCor <- maxCorPPMCD(x, y, alpha=alpha, nIterations=nIterations, 
                          nAlternate=nAlternate, nGrid=nGrid, tol=tol, 
                          seed=seed)
  }
  ## return maximum correlation
  maxCor$call <- matchedCall
  class(maxCor) <- "maxCor"
  maxCor
}


#' Robust CCA via alternating series of grid searches based on the bivariate MCD
#' 
#' Perform robust canoncial correlation analysis via projection pursuit based 
#' on alternating series of grid searches in two-dimensional subspaces of each 
#' data set.  The correlations in those subspaces are computed with the 
#' bivariate MCD-estimator.
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
#' The order of the variables in a series of grid searches for each of the data 
#' sets is determined by the average absolute correlations with the variables 
#' of the respective other data set.  This requires to compute the full 
#' \eqn{(p \times q)}{(p x q)} matrix of absolute correlations, where \eqn{p} 
#' denotes the number of variables of \code{x} and \eqn{q} the number of 
#' variables of \code{y}.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
#' @param alpha  a numeric value controlling the percentage of observations 
#' over which the determinant is minimized.  Allowed values are between 0.5 
#' and 1 and the default is 0.5.
#' @param nIterations  an integer giving the maximum number of iterations.
#' @param nAlternate  an integer giving the maximum number of alternate series 
#' of grid searches in each iteration.
#' @param nGrid  an integer giving the number of equally spaced grid points on 
#' the unit circle to use in each grid search.
#' @param tol  a small positive numeric value to be used for determining 
#' convergence.
#' @param fallback  logical; the data are first robustly standardized via 
#' median and MAD.  This indicates whether standardization via mean and 
#' standard deviation should be performed as a fallback mode for variables 
#' whose MAD is zero (e.g., for dummy variables).
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
#' @note Unlike \code{\link{ccaGrid}}, \code{ccaMCD} is a pure \R implementation 
#' of the alternate grid algorithm.  For this reason, but mostly due to the 
#' computational cost of the MCD algorithm, \code{ccaMCD} is much slower than 
#' \code{\link{ccaGrid}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaGrid}}, \code{\link{maxCorMCD}}, \code{\link{corMCD}}
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
#' ## robust CCA
#' ccaMCD(x, y, alpha = 0.75)
#' 
#' @keywords multivariate robust
#' 
#' @export

ccaMCD <- function(x, y, k = 1, alpha = 0.5, nIterations = 10, 
                   nAlternate = 10, nGrid = 25, tol = 1e-06, 
                   standardize = TRUE, fallback = FALSE) {
  ## initializations
  matchedCall <- match.call()
  x <- as.matrix(x)
  y <- as.matrix(y)
  p <- ncol(x)
  q <- ncol(y)
  k <- as.integer(k)
  standardize <- isTRUE(standardize)
  fallback <- isTRUE(fallback)
  ## check control arguments for PP algorithm
  nIterations <- as.integer(nIterations)
  nAlternate <- as.integer(nAlternate)
  nGrid <- as.integer(nGrid)
  tol <- as.numeric(tol)
  ## check control arguments for correlation based on bivariate MCD-estimator
  alpha <- as.numeric(alpha)
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) runif(1)
  seed <- .Random.seed
  on.exit(set.seed(seed))  # reset seed of the random number generator
  ## initialize results
  r <- numeric(k)
  A <- matrix(0, p, k)
  B <- matrix(0, q, k)
  ## robustly standardize data if requested
  if(standardize) {
    xs <- robStandardize(x, fallback=fallback)
    ys <- robStandardize(y, fallback=fallback)
    scaleX <- attr(xs, "scale")
    scaleY <- attr(ys, "scale")
    # compute first canonical correlation variables with standardized data
    # and transform canonical vectors back to original scale
    tmp <- maxCorPPMCD(xs, ys, alpha=alpha, nIterations=nIterations, 
                       nAlternate=nAlternate, nGrid=nGrid, tol=tol, 
                       seed=seed)
    a <- backtransform(tmp$a, scaleX)
    b <- backtransform(tmp$b, scaleY)
  } else {
    # compute first canonical correlation variables with original data
    tmp <- maxCorPPMCD(x, y, alpha=alpha, nIterations=nIterations, 
                       nAlternate=nAlternate, nGrid=nGrid, tol=tol, 
                       seed=seed)
    a <- tmp$a
    b <- tmp$b
  }
  ## store first canonical variables
  r[1] <- tmp$cor
  A[, 1] <- a
  B[, 1] <- b
  ## compute remaining canonical correlations
  if(k > 1) {
    xl <- x
    yl <- y
    for(l in 2:k) {
      # perform Householder transformation
      Pl <- householder(a)
      Ql <- householder(b)
      xl <- (xl %*% Pl)[, -1, drop=FALSE]  # reduced x data
      yl <- (yl %*% Ql)[, -1, drop=FALSE]  # reduced y data
      # standardize transformed data if requested
      if(standardize) {
        xs <- robStandardize(xl, fallback=fallback)
        ys <- robStandardize(yl, fallback=fallback)
        scaleX <- attr(xs, "scale")
        scaleY <- attr(ys, "scale")
        # compute the canonical correlation and canonical vectors for the
        # standardized reduced data sets and transform the canonical
        # vectors back to the original scale
        tmp <- maxCorPPMCD(xs, ys, alpha=alpha, nIterations=nIterations, 
                           nAlternate=nAlternate, nGrid=nGrid, tol=tol, 
                           seed=seed)
        a <- backtransform(tmp$a, scaleX)
        b <- backtransform(tmp$b, scaleY)
      } else {
        # compute canonical correlation and canonical vectors for reduced data
        tmp <- maxCorPPMCD(xl, yl, alpha=alpha, nIterations=nIterations, 
                           nAlternate=nAlternate, nGrid=nGrid, tol=tol, 
                           seed=seed)
        a <- tmp$a
        b <- tmp$b
      }
      # extract canonical correlation
      r[l] <- tmp$cor
      # transform canonical vectors back to original space
      if(l == 2) {
        P <- Pl
        Q <- Ql
      } else {
        # expand current Householder matrix for x and premultiply with
        # product of previous ones
        Pl <- rbind(matrix(0, l-2, p-l+2), Pl)
        Pl <- cbind(diag(nrow=p, ncol=l-2), Pl)
        P <- P %*% Pl
        # expand current Householder matrix for y and premultiply with
        # product of previous ones
        Ql <- rbind(matrix(0, l-2, q-l+2), Ql)
        Ql <- cbind(diag(nrow=q, ncol=l-2), Ql)
        Q <- Q %*% Ql
      }
      # expand canonical vectors and premultiply with product of
      # corresponding Householder matrices
      A[, l] <- P %*% c(rep.int(0, l-1), a)
      B[, l] <- Q %*% c(rep.int(0, l-1), b)
    }
  }
  ## return results
  cca <- list(cor=r, A=A, B=B, call=matchedCall)
  class(cca) <- "cca"
  cca
}


## workhorse function for maximum correlation via alternate grid searches
maxCorPPMCD <- function(x, y, alpha = 0.5, nIterations = 10, 
                        nAlternate = 10, nGrid = 25, 
                        tol = 1e-06, seed = NULL) {
  # initializations
  p <- ncol(x)
  q <- ncol(y)
  # perform alternate grid searches if both data sets are multivariate
  # if one data set is univariate, alternate grid searches are not necessary
  if(p == 1 && q == 1) {
    # both data sets are univariate
    a <- b <- 1
    maxCor <- abs(corMCD(x, y, alpha=alpha, seed=seed))  # compute correlation
  } else {
    previousMaxCor <- -Inf
    if(p > 1 && q == 1) {
      # x is multivariate, y is univariate
      tmp <- findOrder(x, y, alpha=alpha, seed=seed)  # column order
      orderX <- tmp$orderX
      maxCor <- tmp$maxCor
      a <- tmp$a
      b <- 1
      # stop if there are two consecutive iterations without improvement
      i <- 0
      convCounter <- 0
      while(i < nIterations && convCounter < 2) {
        previousMaxCor <- maxCor
        grid <- getGrid(nGrid, i+1)  # define vector of grid points
        tmp <- gridSearch(x, orderX, y, alpha, seed, grid, maxCor, a)
        maxCor <- tmp$maxCor
        a <- tmp$a
        i <- i + 1
        convCounter <- setCounter(convCounter, maxCor, previousMaxCor, tol)
      }
    } else if(p == 1 && q > 1) {
      # x is univariate, y is multivariate
      tmp <- findOrder(y, x, alpha=alpha, seed=seed)  # column order
      orderY <- tmp$orderX
      maxCor <- tmp$maxCor
      b <- tmp$a
      a <- 1
      # stop if there are two consecutive iterations without improvement
      i <- 0
      convCounter <- 0
      while(i < nIterations && convCounter < 2) {
        previousMaxCor <- maxCor
        grid <- getGrid(nGrid, i+1)  # define vector of grid points
        tmp <- gridSearch(y, orderY, x, alpha, seed, grid, maxCor, b)
        maxCor <- tmp$maxCor
        b <- tmp$a
        i <- i + 1
        convCounter <- setCounter(convCounter, maxCor, previousMaxCor, tol)
      }
    } else if(p > 1 && q > 1) {
      # both data sets are multivariate
      # compute correlations between variables in x and y
      tmp <- findOrder(x, y, alpha=alpha, seed=seed)  # column orders
      orderX <- tmp$orderX
      orderY <- tmp$orderY
      maxCor <- tmp$maxCor
      a <- tmp$a
      b <- tmp$b
      # perform alternate grid searches
      if(tmp$startWithX) {
        # start with grid search for x
        # stop if there are two consecutive iterations without improvement
        i <- 0
        convCounter <- 0
        while(i < nIterations && convCounter < 2) {
          previousMaxCor <- maxCor
          grid <- getGrid(nGrid, i+1)  # define vector of grid points
          j <- 0
          altMaxCor <- -Inf
          while(j < nAlternate && (maxCor - altMaxCor) > tol) {
            altMaxCor <- maxCor
            # maximize correlation functional over a keeping b fixed
            tmp <- gridSearch(x, orderX, y %*% b, alpha, seed, grid, maxCor, a)
            maxCor <- tmp$maxCor
            a <- tmp$a
            # maximize correlation functional over b keeping a fixed
            tmp <- gridSearch(y, orderY, x %*% a, alpha, seed, grid, maxCor, b)
            maxCor <- tmp$maxCor
            b <- tmp$a
            # update iterator
            j <- j + 1
          }
          # update iterators
          i <- i + 1
          convCounter <- setCounter(convCounter, maxCor, previousMaxCor, tol)
        }
      } else {
        # start with grid search for y
        # stop if there are two consecutive iterations without improvement
        i <- 0
        convCounter <- 0
        while(i < nIterations && convCounter < 2) {
          previousMaxCor <- maxCor
          grid <- getGrid(nGrid, i+1)  # define vector of grid points
          j <- 0
          altMaxCor <- -Inf
          while(j < nAlternate && (maxCor - altMaxCor) > tol) {
            altMaxCor <- maxCor
            # maximize correlation functional over b keeping a fixed
            tmp <- gridSearch(y, orderY, x %*% a, alpha, seed, grid, maxCor, b)
            maxCor <- tmp$maxCor
            b <- tmp$a
            # maximize correlation functional over a keeping b fixed
            tmp <- gridSearch(x, orderX, y %*% b, alpha, seed, grid, maxCor, a)
            maxCor <- tmp$maxCor
            a <- tmp$a
            # update iterator
            j <- j + 1
          }
          # update iterators
          i <- i + 1
          convCounter <- setCounter(convCounter, maxCor, previousMaxCor, tol)
        }
      }
    }
  }
  # ensure that norm of weighting vectors is 1
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  # check direction
  r <- corMCD(x %*% a, y %*% b)
  if(r < 0) {
    if(p > 1 && q == 1) a <- -a
    else b <- -b
  }
  # return maximum correlation
  list(cor=maxCor, a=a, b=b)
}


## utility functions

# (robustly) standardize the data
standardize <- function(x, robust = TRUE) {
  if(robust) {
    # with median and MAD
    tmp <- apply(x, 2, function(v) unlist(fastMAD(v)))
    center <- tmp[1,]  # column medians
    x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
    scale <- tmp[2,]  # column MADs
    x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
  } else {
    # with mean and standard deviation
    center <- colMeans(x)  # compute column means (faster than apply)
    x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
    f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
    scale <- apply(x, 2, f)  # compute column scales with zero means
    x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
  }
  # add attributes and return standardized data
  attr(x, "center") <- center
  attr(x, "scale") <- scale
  x
}

# robustly standardize the data with fallback option
robStandardize <- function(x, fallback = FALSE) {
  # robustly standardize data
  xs <- standardize(x, robust=TRUE)
  # if requested, check if some variables have a robust scale of zero and 
  # use standardization with mean/sd as fallback mode
  if(fallback) {
    scale <- attr(xs, "scale")
    if(ncol(x) == 1) {
      if(scale == 0) xs <- standardize(x)
    } else {
      isZero <- which(scale == 0)
      if(length(isZero) > 0) {
        # standardize with mean and standard deviation
        center <- attr(xs, "center")
        xcs <- standardize(x[, isZero])
        center[isZero] <- attr(xcs, "center")
        scale[isZero] <- attr(xcs, "scale")
        xs[, isZero] <- xcs
        attr(xs, "center") <- center
        attr(xs, "scale") <- scale
      }
    }
  }
  # return standardized data
  xs
}

# get the order in which to update the elements of a weighting vector
findOrder <- function(x, y, alpha = 0.5, seed = NULL) {
  # initializations
  p <- ncol(x)
  q <- ncol(y)
  # get order
  if(q == 1) {
    # compute columnwise absolute correlations of x with y
    corY <- abs(apply(x, 2, corMCD, y, alpha=alpha, seed=seed))
    # order columns of x according to absolute correlations with y
    orderX <- order(corY)
    # find maximum correlation and set weight of corresponding variable to 1
    first <- orderX[1]
    a <- rep.int(0, p)
    a[first] <- 1
    # return order, maximum correlation and weighting vector as list
    list(orderX=orderX, maxCor=corY[first], a=a)
  } else {
    # compute complete matrix of columnwise absolute correlations
    corMat <- abs(apply(y, 2, function(y) {
      apply(x, 2, corMCD, y, alpha=alpha, seed=seed)
    }))
    # compute average absolute correlations
    avgCorX <- rowMeans(corMat)
    avgCorY <- colMeans(corMat)
    # determine the order of the variables and the starting values for the
    # weighting vectors from average correlations
    orderX <- order(avgCorX)
    orderY <- order(avgCorY)
    # set weights of corresponding variables to 1 and determine with which
    # data set to start in alternate grid searches
    firstX <- orderX[1]
    a <- rep.int(0, p)
    a[firstX] <- 1
    firstY <- orderY[1]
    b <- rep.int(0, q)
    b[firstY] <- 1
    # return orders, maximum correlation and weighting vectors as list
    list(orderX=orderX, orderY=orderY, maxCor=corMat[firstX, firstY], 
         a=a, b=b, startWithX=avgCorX[firstX] >= avgCorY[firstY])
  }
}

# get equispaced grid of angles in a plane
getGrid <- function(nGrid = 25, i = 1) {
  seq(-pi/2^i, by=pi/(nGrid*2^(i-1)), length.out=nGrid)
}

# grid search to update one weighting vector
gridSearch <- function(x, orderX, y, alpha = 0.5, seed = NULL, grid, maxCor, a) {
  ## initializations
  p <- ncol(x)
  ## perform grid searches for each canonical basis vector
  for(j in seq_len(p)) {
    # define current canonical basis vector according to order of columns
    ej <- rep.int(0, p)
    ej[orderX[j]] <- 1
    # perform grid search for the current canonical basis vector
    corY <- abs(sapply(grid, function(angle) {
      currentA <- cos(angle) * a + sin(angle) * ej;
      corMCD(x %*% currentA, y, alpha=alpha, seed=seed)
    }))
    # find grid point that maximizes correlation functional and keep
    # maximum correlation of current grid search
    whichMax <- which.max(corY)
    currentMaxCor <- corY[whichMax]
    # update maximum correlation and weighting vector
    # if 0 degree angle is not part of the grid, the maximum correlation
    # of the current grid search may be smaller than the previous maximum
    if(currentMaxCor > maxCor) {
      maxCor <- currentMaxCor
      optAngle <- grid[whichMax]
      a <- cos(optAngle) * a + sin(optAngle) * ej
    }
  }
  # return maximum correlation and weighting vector
  list(maxCor=maxCor, a=a)
}

# set counter of consecutive iterations without improvement
setCounter <- function(counter, current, previous, tol = 1e-06) {
  if((current - previous) > tol) 0 else counter + 1
}


# transform canonical vectors back to the original scale
backtransform <- function(a, scale) {
  a <- a / scale           # divide by scale of corresponding variable
  a <- a / sqrt(sum(a^2))  # divide by norm
}

# compute rotation matrix for Householder transformation
householder <- function(a) {
  p <- length(a)
  e1 <- c(1, rep.int(0, p-1))  # first basis vector
  n1 <- e1 - a                 # normal vector
  n1 <- n1 / sqrt(sum(n1^2))    # unit normal vector
  diag(p) - 2 * tcrossprod(n1)
}
