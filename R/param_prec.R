#' Precision matrix parametrization helper functions.
#' @rdname prec
#' @param theta numeric vector of length `m`.
#' @param p numeric giving the dimention of Q.
#' If missing, `p = (1+sqrt(1+8*length(theta)))`
#' and Q is assumed to be dense.
#' @param ilowerL numeric vector as index
#' to (lower) L to be filled with `theta`.
#' Default is missing and Q is assumed to be dense.
#' @return matrix as the Cholesky factor of a
#' precision matrix as the inverse of a correlation
#' @details
#' The precision matrix definition consider
#' `m` parameters for the lower part of L.
#' If Q is dense, then `m = p(p-1)/2`, else
#' `m = length(ilowerL)`.
#' Then the precision is defined as
#' \eqn{Q(\theta) = L(\theta)L(\theta)^T}
#' @return a matrix whose elements at the lower
#' triangle are the filled in elements of the
#' Cholesky decomposition of a precision matrix
#' and diagonal elements as `1:p`.
#' @export
#' @examples
#' theta1 <- c(1, -1, -2)
#' Lprec(theta1)
#' theta2 <- c(0.5, -0.5, -1, -1)
#' Lprec(theta2, 4, c(2,4,7,12))
Lprec <- function(theta, p, ilowerL) {
  if(missing(p)) {
    m <- length(theta)
    p <- (1 + sqrt(1+8*m))/2
    stopifnot((p==floor(p)) & (p==ceiling(p)))
  }
  stopifnot(p>0)
  L <- diag(x=1:p, nrow = p, ncol = p)
  if(missing(ilowerL)) {
    L[lower.tri(L)] <- theta
  } else {
    L[ilowerL] <- theta
    G <- diag(p)
    G[ilowerL] <- -1
    G <- t(G)
    G[ilowerL] <- -1
    diag(G) <- 1 - colSums(G)
    lG <- t(chol(G))
    ifill <- which(is.zero(G) & (!is.zero(lG)))
    if(length(ifill)>0) {
      L <- fillLprec(L, ifill)
    }
  }
  return(L)
}
#' @describeIn prec
#' Function to fill-in a Cholesky matrix
#' @param L matrix as the lower triangle
#' containing the Cholesky decomposition of
#' a precision matrix
#' @param lfi indicator of fill-in elements
fillLprec <- function(L, lfi) {
  L <- as.matrix(L)
  ii <- row(L)[lfi]
  jj <- col(L)[lfi]
  for(v in 1:length(ii)) {
    i <- ii[v]
    j <- jj[v]
    if(j==1) {
      warning("j = 1!\n")
      L[i,1] <- 0.0
    } else {
      stopifnot((i>1) & (j>1)) ## L_{11} not allowed
      stopifnot(j>1) ## j=1 is not allowed
      k <- 1:(j-1)
      L[i, j] <- -sum(L[i, k] * L[j, k]) / L[j, j]
    }
  }
  return(L)
}

