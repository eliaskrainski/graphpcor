#' Functions used by/to basecor
#' @name basecor-utils
NULL
#> NULL

#' @describeIn basecor-utils
#' Build the (lower) Cholesky for a correlation matrix
#' @inheritParams basecor
#' @param theta numeric parameter vector.
#' @returns matrix with lower triangle as the Cholesky factor
#' of a correlation matrix, and atttibutes:
#' `parametrization` ("cpc" or "sap"), "theta"
#' the parameter vector under the parametrization,
#' "iLtheta" the index in the lower triangle matrix for theta,
#' "logDeterminant" the log determinant for the correlation matrix.
#' @export
cholcor <- function(
    theta,
    p,
    iLtheta,
    parametrization = "cpc") {
  parametrization <- match.arg(
    arg = tolower(parametrization),
    choices = c("cpc", "sap")
  )
  stopifnot((m <- length(theta))>0)
  if(missing(p)) {
    p <- (1 + sqrt(1+8*m))/2
  }
  stopifnot(floor(p)==ceiling(p))
  stopifnot(p>1)
  ith0 <- which(lower.tri(
    diag(x = rep(1, p), nrow = p, ncol = p)))
  if(missing(iLtheta)) {
    iLtheta <- ith0
  } else {
    stopifnot(all(iLtheta %in% ith0))
  }
  stopifnot(length(theta)==length(iLtheta))
  B <- A <- diag(p)
  if(parametrization == 'cpc') {
    A[iLtheta] <- tanh(theta)
    B[ith0] <- sqrt(1-A[ith0]^2)
  } else {
    theta <- pi/(1+exp(-theta))
    A[iLtheta] <- cos(theta)
    B[ith0] <- 1.0
    B[iLtheta] <- sin(theta)
  }
  if(p>2) {
    for(j in 2:(p-1)) {
      B[, j] <- B[, j] * B[, j-1]
    }
  }
  L <- A * cbind(1, B[, 1:(p-1)])
  attr(L, "parametrization") <- parametrization
  attr(L, "theta") <- theta
  attr(L, "iLtheta") <- iLtheta
  attr(L, "logDeterminant") <- sum(diag(L))*2.0
  return(L)
}
#' Map a correlation matrix to theta.
#' @describeIn basecor-utils
#' This function takes a correlation matrix and return
#' the parameter vector that generates it with the
#' canonical partial correlation parametrization.
#' @param corr matrix as a correlation matrix
corr2CPC_theta <- function(corr) {
  corr <- as.matrix(corr)
  stopifnot(nrow(corr) == (p <- ncol(corr)))
  stopifnot(p>1)
  stopifnot(all.equal(
    new("numeric", diag(corr)),
    rep(1.0, p)))
  stopifnot(all.equal(corr,t(corr)))
  L <- t(chol(corr))
  B <- A <- diag(p)
  A[, 1] <- L[, 1]
  if(p>2) {
    B[2:p, 1] <- sqrt(1 - (A[2:p, 1]^2))
    for(j in 2:(p-1)) {
      A[j:p, j] <- L[j:p, j]/B[j:p, j-1]
      sj <- pmax(0, 1 - (A[j:p, j]^2))
      B[j:p, j] <- sqrt(sj) * B[j:p, j-1]
    }
  }
  atanh(A[which(lower.tri(L))])
}
#' Draw samples from a `basecor`.
#' @describeIn basecor-utils
#' Sample from the model parameters and map it to the correlation matrix.
#' @param x the correlation model
#' @param size the number of samples
#' @param lambda the penalization parameter
#' @export
sample.basecor <- function(x, size, lambda) {
  stopifnot((m <- length(x$theta))>0)
  stopifnot(lambda>0)
  stopifnot(size>0)
  p <- ncol(x$base)
  r <- rexp(size, lambda)
  theta <- t(sapply(1:size, function(i) {
    z <- rnorm(m)
    z <- z / sqrt(sum(z^2))
  })) * r
  H <- hessian(x)
  sHi <- graphpcor:::dspd(H)$sqrtInv
  theta <- sweep(theta %*% sHi, 2, x$theta)
  out <- sapply(1:size, function(i){
    tcrossprod(cholcor(theta[i, ], p))
  })
  dim(out) <- c(p, p, size)
  return(out)
}
