#' Precision matrix parametrization helper functions.
#' @rdname prec
#' @param theta numeric vector of length `m`
#' with the parameter
#' @param p integer giving the dimention of Q.
#' If `p` and `ilowerL` are missing, then
#' Q is assumed to be dense and
#' `p = (1+sqrt(1+8*length(theta)))`.
#' @param ilowerL integer vector as index
#' to (lower) L to be filled with `theta`.
#' @return matrix as the Cholesky factor of a
#' precision matrix as the inverse of a correlation
#' @details
#' The precision matrix definition consider
#' `m` parameters for the lower part of L.
#' If Q is dense, then `m = p(p-1)/2`, else
#' `m = length(ilowerL)`.
#' The precision is defined as
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
Lprec <- function(theta, p, ilowerL, d0) {
  if(missing(p)) {
    if(!missing(ilowerL)) {
      stop("Please provide `p`!")
    }
    m <- length(theta)
    p <- (1 + sqrt(1+8*m))/2
    stopifnot((p==floor(p)) & (p==ceiling(p)))
  }
  stopifnot(p>0)
  if(missing(d0)) {
    d0 <- 1:p
  }
  stopifnot(length(d0)==p)
  L <- diag(x=d0, nrow = p, ncol = p)
  if(missing(ilowerL)) {
    L[lower.tri(L)] <- theta
  } else {
    L[ilowerL] <- theta
    L <- fillLprec(L, ilowerL = ilowerL)
  }
  return(L)
}
#' @describeIn prec
#' Function to fill-in a Cholesky matrix
#' @param L matrix as the lower triangle
#' containing the Cholesky decomposition of
#' a precision matrix
#' @param lfi integer vector used as indicator of the
#' position in the lower matrix where are the
#' fill-in elements. Must be col then row ordered.
fillLprec <- function(L, lfi, ilowerL) {
  L <- as.matrix(L)
  p <- nrow(L)
  if(missing(lfi)) {
    G <- matrix(0, p, p)
    G[ilowerL] <- -1
    G <- t(G)
    G[ilowerL] <- -1
    diag(G) <- 1 - colSums(G)
    lG <- t(chol(G))
    lfi <- which(is.zero(G) & (!is.zero(lG)))
  }
  if(length(lfi)>0) {
    if(length(lfi)>1)
      stopifnot(all(diff(lfi)>0))
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
  }
  return(L)
}
#' @describeIn prec
#' Build a correlation matrix from the
#' precision parametrization
theta2correl <- function(theta, p, ilowerL, d0,
                         transf = c("natural", "cpc")) {
  transf <- match.arg(transf, c("natural", "cpc"))
  if(transf=="natural") {
    L <- Lprec(theta, p=p, ilowerL=ilowerL, d0=d0)
    V <- chol2inv(t(L))
    return(cov2cor(V))
  } else {
    crossprod(Lcorrel(theta, p=p, ilowerL=ilowerL))
  }
}
#' @describeIn correl
#' Cholesky of a correlation matrix using the
#'  [correlation-matrix-inverse-transform](https://mc-stan.org/docs/reference-manual/transforms.html)
#' where tanh(\eqn{\theta_j}) is the canonical
#' partial correlation - CPC.
#' @param ilowerL integer vector as index
#' to (lower) L to be filled with `theta`.
#' Assume dense if missing.
#' @export
Lcorrel <- function(theta, p, ilowerL) {
  if(missing(p)) {
    if(!missing(ilowerL)) {
      stop("Please provide `p`!")
    }
    m <- length(theta)
    p <- (1 + sqrt(1+8*m))/2
    stopifnot((p==floor(p)) & (p==ceiling(p)))
  }
  stopifnot(p>0)
  z <- diag(x=rep(1,p), nrow = p, ncol = p)
  if(missing(ilowerL)) {
    z[lower.tri(z)] <- tanh(theta)
  } else {
    z[ilowerL] <- tanh(theta)
  }
  z <- fillLprec(z, ilowerL = ilowerL)
  w <- z <- t(z)
  psz <- rep(1, p)
  for(i in 2:p) {
    psz <- psz * sqrt(1-z[i-1, ]^2)
    w[i,i:p] <- z[i,i:p] * psz[i:p]
  }
  return(w)
}
