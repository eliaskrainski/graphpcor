#' Build the Cholesky factor, L, of a precision matrix.
#' @param theta numeric vector of length `m`
#' with the parameters.
#' @param p integer giving the dimension of
#' the precision matrix.
#' @param itheta integer vector as index
#' to (the lower part of) L
#' to be filled with `theta`.
#' Length of itheta must be equal 'm'.
#' @param d0 elements at the diagonal of L.
#' By default uses `p:1`.
#' @return matrix as the (lower triangle)
#' Cholesky factor of a precision matrix.
#' @details
#' The precision is defined as
#' \eqn{Q(\theta) = L(\theta)L(\theta)^T}
#' where \eqn{L(\theta)} the Cholesky of \eqn{Q(\theta)}
#' is filled using the following steps
#' 1. Build \eqn{L} starting with filling its diagonal
#'    with 'd0'.
#' 2. Place \eqn{theta} at the positions 'itheta'
#' 3. Compute the fill-in elements
#' @seealso [correl]
#' @export
#' @examples
#' Lprec(c(1,-1,-2), 3, c(2,3,6))
#' Lprec(c(2,-1,1,-0.5), 4, c(2,3,8,12))
Lprec <- function(theta, p, itheta, d0) {
  stopifnot(p>1)
  if(missing(d0)) {
    d0 <- p:1
  }
  stopifnot(length(d0)==p)
  L <- diag(x=d0, nrow = p, ncol = p)
  if(missing(itheta)) {
    stopifnot(length(theta)==(p*(p-1)/2))
    L[lower.tri(L)] <- theta
  } else {
    L[itheta] <- theta
    L <- fillLprec(L)
  }
  return(L)
}
#' @describeIn Lprec
#' Function to fill-in a Cholesky matrix
#' @param L matrix as the lower triangle
#' containing the Cholesky decomposition of
#' a precision matrix
#' @param lfi integer vector used as indicator of the
#' position in the lower matrix where are the
#' fill-in elements. Must be col then row ordered.
fillLprec <- function(L, lfi) {
  L <- as.matrix(L)
  p <- nrow(L)
  if(missing(lfi)) {
    i0 <- is.zero(L)
    G <- i0 - 1.0
    G <- G + t(G)
    diag(G) <- 1 - colSums(G)
    lG <- t(chol(G))
    lfi <- which(i0 & (!is.zero(lG)))
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
