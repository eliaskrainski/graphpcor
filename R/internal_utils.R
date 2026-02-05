#' Internal functions.
#' @name internal-utils
NULL
#> NULL

#' @describeIn internal-utils
#' Function to deal with `p` and `itheta`
#' @inheritParams basepcor
p_itheta_fncheck <- function(p, itheta) {
  if(missing(itheta)) {
    if(missing(p))
      stop("Please provide 'p' or 'itheta'!")
    p <- as.integer(p[1])
    itheta <- which(lower.tri(diag(
      x = rep(1, p), nrow = p, ncol = p)))
  } else {
    if(inherits(itheta, "graphpcor")) {
      Q1 <- Laplacian(itheta)
      p <- ncol(Q1)
      itheta <- which(lower.tri(Q1) & (Q1 != 0.0))
    } else {
      if(missing(p))
        stop("Please provide 'p'!")
    }
  }
  il0 <- which(lower.tri(matrix(0, p, p)))
  stopifnot(all(itheta %in% il0))
  attr(itheta, 'p') <- as.integer(p)
  return(itheta)
}
#' @describeIn internal-utils
#' Function to deal with `m` and `iparams`
#' @param m integer to specify the number of parameters
m_iparams_fncheck <- function(m, iparams) {
  if(missing(iparams) || is.null(iparams)) {
    if(missing(m)) {
      stop("Missing 'm' and 'iparams'!")
    }
    iparams <- 1:m
  }
  stopifnot(length(iparams) == m)
  stopifnot(all(iparams %in% 1:m))
  ## next text allow c(1,1,2,2,1) but not c(2,2,1,1,2)
  stopifnot(all(diff(sort(unique(iparams)))==1))
  attr(iparams, "m") <- m
  return(iparams)
}
#' @describeIn internal-utils
#' Compute the KLD between two multivariate Gaussian
#' distributions, assuming equal mean vector
#' @param C1 is a correlation matrix.
#' @param C0 is a correlation matrix of the base model.
#' @param L1 is the Cholesky of `C1`.
#' @param L0 is the Cholesky of `C0`.
#' @details
#' By assuming equal mean vector we have
#'  \deqn{KLD = 0.5( tr(C0^{-1}C1) -p - log(|C1|) + log(|C0|) )}
KLD10 <- function(C1, C0, L1, L0) {
  ### input: C1, C0 or, alternatively, L1, L0 (upper triangles)
  ### output: KLD
  if(missing(L1)) {
    if(missing(C1)) {
      stop("Please provide either 'C1' or 'L1'!")
    }
    L1 <- chol(C1)
  }
  if(missing(C1)) {
    C1 <- crossprod(L1)
  }
  p <- nrow(L1)
  hld1 <- sum(log(diag(L1)))
  if(missing(C0)) {
    if(missing(L0)) {
      warning("Missing C0,L0: using 'I'!")
      L0 <- diag(x = rep(1, p), nrow = p, ncol = p)
    }
  } else {
    if(missing(L0)) {
      L0 <- chol(C0)
    }
  }
  hld0 <- sum(log(diag(L0)))
  tr <- sum(diag(chol2inv(L0) %*% C1))
  return(0.5*(tr -p) + hld0 - hld1)
}
#' @describeIn internal-utils
#' Function (internal) to decompose a positive definite matrix,
#' and compute useful elements out of that.
#' @param x matrix.
#' @param decomposition character to inform
#' which decomposition is to be applied to the
#' hessian. The options are "eigen", "svd" and "chol".
#' Default is "svd".
#' @returns `dspd` returns a list with the decomposition elements,
#'  "logDeterminant" (of the original matrix),
#' "sqrt" (its 'square root') and
#' "sqrtInv" (its inverse 'square root').
dspd <- function(x, decomposition = "svd") {
  decomposition <- match.arg(
    tolower(decomposition),
    c("svd", "eigen", "chol"))
  stopifnot(nrow(x)==ncol(x))
  p <- ncol(x)
  ## next bit follows mvtnorm:::rmvnorm()
  t0 <- sqrt(.Machine$double.eps)
  if(decomposition == "eigen") {
    xd <- eigen(x)
    tol1 <- t0 * abs(xd$values[1])
    if(!all(xd$values >= tol1))
      warning("'x' is numerically not positive semidefinite")
    xd$logDeterminant <- sum(log(xd$values))
    s <- sqrt(pmax(xd$values, 0.0))
    xd$sqrt <- t(xd$vectors %*% (t(xd$vectors) * s))
    xd$sqrtInv <- t(xd$vectors %*% (t(xd$vectors) / s))
  }
  if(decomposition == "svd") {
    xd <- svd(x)
    tol1 <- t0 * abs(xd$d[1])
    if(any(xd$d<tol1))
      warning("'x' is numerically not positive semidefinite")
    xd$logDeterminant <- sum(log(xd$d))
    s <- sqrt(pmax(xd$d, 0.0))
    xd$sqrt <- t(xd$v %*% (t(xd$u) * s))
    xd$sqrtInv <- t(xd$v %*% (t(xd$u) / s))
  }
  if(decomposition == "chol") {
    xd <- chol(x, pivot = TRUE)
    if(any(diag(xd)<t0))
      warning("'x' is numerically not positive semidefinite")
    tol1 <- pmin(t0 * 10, diag(xd))
    xd$logDeterminant <- sum(diag(xd))*2.0
    x$sqrt <- matrix(xd[, order(attr(xd, "pivot")), ], p)
    xn <- chol2inv(chol(x))
    xn.5 <- chol(xn, pivot = TRUE)
    xd$sqrtInv <- matrix(xn.5[, order(attr(xn.5, "pivot"))], p)
  }
  return(xd)
}
