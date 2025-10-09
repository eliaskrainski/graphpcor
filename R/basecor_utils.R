#' @describeIn basecor Cholesky (lower triangular) matrix from theta.
#' @returns matrix with lower triangle as the Cholesky factor
#' of a correlation matrix if parametrization is
#' "cpc" or "sap" and of a precision matrix if
#' parametrization is 'itp' (with 'd0' as the diagonal elements).
theta2L <- function(
    theta,
    p,
    parametrization = "cpc",
    itheta,
    d0) {
  parametrization <- match.arg(
    arg = tolower(parametrization),
    choices = c("cpc", "sap", "itp")
  )
  stopifnot((m <- length(theta))>0)
  if(missing(p)) {
    p <- (1 + sqrt(1+8*m))/2
  }
  stopifnot(floor(p)==ceiling(p))
  stopifnot(p>1)
  ith0 <- which(lower.tri(
    diag(x = rep(1, p), nrow = p, ncol = p)))
  if(missing(itheta)) {
    itheta <- ith0
  } else {
    stopifnot(all(itheta %in% ith0))
  }
  stopifnot(length(theta)==length(itheta))
  if(parametrization == "itp") {
    if(missing(d0)) {
      warning("Using 'd0 = p:1'!")
      d0 <- p:1
    }
    L <- diag(x = d0, nrow = p, ncol = p)
    L[itheta] <- theta
    L <- fillLprec(L)
  } else {
    B <- A <- diag(p)
    if(parametrization == 'cpc') {
      A[itheta] <- tanh(theta)
      B[ith0] <- sqrt(1-A[ith0]^2)
    } else {
      theta <- pi/(1+exp(-theta))
      A[itheta] <- cos(theta)
      B[ith0] <- 1.0
      B[itheta] <- sin(theta)
    }
    if(p>2) {
      for(j in 2:(p-1)) {
        B[, j] <- B[, j] * B[, j-1]
      }
    }
    L <- A * cbind(1, B[, 1:(p-1)])
  }
  attr(L, "parametrization") <- parametrization
  attr(L, "theta") <- theta
  attr(L, "itheta") <- itheta
  if(parametrization == 'itp') {
    attr(L, "d0") <- d0
  }
  return(L)
}

#' Internal functions to basecor
#' @name basecor-utils
NULL
#> NULL

#' @describeIn basecor-utils
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
#' @describeIn basecor-utils
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
  ### input: C1, C0 (alternatively L1, L0)
  ### output: KLD
  if(missing(L1)) {
    if(missing(C1)) {
      stop("Please provide either 'C1' or 'L1'!")
    }
    L1 <- chol(C1)
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
#' @describeIn basecor-utils
#' Evaluate the hessian of the KLD for a
#' correlation model around a base model.
#' @inheritParams basecor
#' @return list containing the hessian,
#' its 'square root', inverse 'square root' along
#' with the decomposition used
#' @importFrom stats cov2cor
#' @importFrom numDeriv hessian
Hcorrel <- function(
    theta,
    p,
    parametrization,
    itheta,
    d0,
    C0,
    decomposition,
    ...) {

  theta2correl <- function(th) {
    L <- theta2L(
      theta = th,
      p = p,
      parametrization = parametrization,
      itheta = itheta,
      d0 = d0
    )
    if("parametrization" == 'itp') {
      return(cov2cor(chol2inv(t(L))))
    } else {
      return(tcrossprod(L))
    }
  }
  if(missing(C0)) {
    C0 <- theta2correl(theta)
  }
  L0 <- chol(C0)
  H <- hessian(
    func = function(x)
      KLD10(theta2correl(x), L0=L0),
    x = theta,
    ...)
  ## next bit follows mvtnorm:::rmvnorm()
  t0 <- sqrt(.Machine$double.eps)
  if(decomposition == "eigen") {
    Hd <- eigen(H)
    if(!all(Hd$values >= (t0 * abs(Hd$values[1]))))
      warning("'H' is numerically not positive semidefinite")
    s <- sqrt(pmax(Hd$values, 0.0))
    h.5 <- t(Hd$vectors %*% (t(Hd$vectors) * s))
    hneg.5 <- t(Hd$vectors %*% (t(Hd$vectors) / s))
  }
  if(decomposition == "svd") {
    Hd <- svd(H)
    if(any(Hd$d<(t0 * abs(Hd$d[1]))))
      warning("'H' is numerically not positive semidefinite")
    s <- sqrt(pmax(Hd$d, 0.0))
    h.5 <- t(Hd$v %*% (t(Hd$u) * s))
    hneg.5 <- t(Hd$v %*% (t(Hd$u) / s))
  }
  if(decomposition == "chol") {
    Hd <- chol(H, pivot = TRUE)
    h.5 <- matrix(Hd[, order(attr(Hd, "pivot")), ], nrow(H))
    hn <- chol2inv(chol(H))
    hn.5 <- chol(hn, pivot = TRUE)
    hneg.5 <- matrix(hn.5[, order(attr(hn.5, "pivot"))], nrow(H))
  }
  stopifnot(all.equal(H, tcrossprod(h.5)))
  attr(H, "base") <- theta
  attr(H, "h.5") <- h.5
  attr(H, "hneg.5") <- hneg.5
  attr(Hd, "decomposition") <- decomposition
  attr(H, "decomposition") <- Hd
  return(H)
}
