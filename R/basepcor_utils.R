#' Functions used by/to basepcor
#' @name basepcor-utils
NULL
#> NULL

#' @describeIn basepcor-utils
#' Function to fill-in a Cholesky matrix
#' @param L matrix as the lower triangle
#' containing the Cholesky decomposition of
#' a initial precision matrix whose non-zeros are
#' only at the position where the lower triangle
#' side of the precision matrix is also non-zero
#' @param lfi integer vector used as indicator of the
#' position in the lower matrix where are the
#' fill-in elements. Must be col then row ordered.
#' @returns lower triangular matrix with the filled-in
#' elements thus `Q0` can be computed.
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
#' @describeIn basepcor-utils
#' Compute the (lower triangle) Cholesky of the initial precision `Q0`.
#' @inheritParams basepcor
#' @param theta numeric, the parameter vector.
#' @returns lower triangular matrix
#' @details The (lower triangle) Cholesky factor
#' of the initial precision for a correlation matrix contains
#' the parameters in the non-zero elements of the lower triangle side
#' of the precision matrix.
#' The filled-in elements are computed from them using [fillLprec()].
Lprec0 <- function(
    theta,
    p,
    iLtheta,
    d0) {
  stopifnot((m <- length(theta))>0)
  iLtheta <- p_iLtheta_fncheck(p, iLtheta)
  p <- attr(iLtheta, "p")
  stopifnot(p>1)
  if(missing(d0)) {
    warning("Using 'd0 = p:1'!")
    d0 <- p:1
  }
  L <- diag(x = d0, nrow = p, ncol = p)
  L[iLtheta] <- theta
  L <- fillLprec(L)
  return(L)
}
#' Map a correlation matrix to theta.
#' @describeIn basepcor-utils
#' This function takes a correlation matrix and return
#' the parameter vector that approximates it under a `graphpcor`.
#' @param corr matrix as a correlation matrix
#' @inheritParams basepcor d0 ...
corr2graphpcor_theta <- function(corr, ..., d0, iLtheta) {
  corr <- as.matrix(corr)
  stopifnot(nrow(corr) == (p <- ncol(corr)))
  stopifnot(p>1)
  stopifnot(all.equal(corr, t(corr)))
  stopifnot(all.equal(
    new("numeric", diag(corr)),
    rep(1.0, p)))
  if(missing(d0) || is.null(d0)) {
    d0 <- p:1
  }
  stopifnot(length(d0)==p)
  if(missing(iLtheta)) {
    have_qgraph <- try(do.call(
      what = "require",
      args = list(package = "qgraph")), silent = TRUE)
    if(inherits(have_qgraph, "try-error")) {
      cat(have_qgraph)
      stop("Please install the 'qgraph' package!")
    }
    ## Fit the sparse partial correlation matrix
    ggmfit <- qgraph::ggmModSelect(S = corr, ...)
    iLtheta <- which((abs(ggmfit$graph)>0) & lower.tri(corr))
    pCorr0 <- diag(p) - ggmfit$graph
    L0 <- t(chol(pCorr0))
    for(i in 1:p)
      L0[i, ] <- (d0[i]/L0[i,i]) * L0[i,]
  } else {
    if(inherits(iLtheta, "graphpcor")) {
      iLtheta <- which(
        (as.matrix(attr(iLtheta, "graph"))>0) &
          lower.tri(diag(p)))
    }
    ## minimize 0.5 * [ trace(C0^{-1} C1) - p - log(|C1|) + log(|C0|) ]
    U0correl <- chol(corr)
    hl0 <- sum(log(diag(U0correl))) ## log(|C0|)/2
    Qbase <- chol2inv(U0correl) ## C0^{-1}
    lQ1 <- diag(d0, p, p) ## working matrix
    ## get initials
    lQ0 <- t(chol(Qbase))
    for(i in 1:p)
      lQ0[i, ] <- (d0[i]/lQ0[i,i]) * lQ0[i,]
    opt <- optim(lQ0[iLtheta], function(x) {
      lQ1[iLtheta] <- x
      C1 <- cov2cor(chol2inv(t(lQ1)))
      r <- sum(diag(Qbase %*% C1))
      return((r-p)/2 + hl0 -sum(diag(chol(C1))))
    })
    L0 <- diag(d0, p, p)
    L0[iLtheta] <- opt$par
    L0 <- fillLprec(L0)
  }
  U0correl <- chol(cov2cor(chol2inv(t(L0))))
  theta <- L0[iLtheta]
  attr(theta, "iLtheta") <- iLtheta
  attr(theta, "L0") <- L0
  attr(theta, "U0correl") <- U0correl
  return(theta)
}
#' Draw samples from a `basepcor`.
#' @describeIn basepcor-utils
#' Sample from the model parameters and map it to the correlation matrix.
#' @param x the correlation model
#' @param size the number of samples
#' @param lambda the penalization parameter
#' @export
sample.basepcor <- function(x, size, lambda) {
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
  lfi <- setdiff(which(abs(x$L0)>0 & lower.tri(x$L0)),
                 x$iLtheta)
  L0 <- diag(x$d0, p, p)
  if(length(lfi)>0) {
    out <- sapply(1:size, function(i){
      L0[x$iLtheta] <- theta[i, ]
      cov2cor(tcrossprod(fillLprec(L0,lfi)))
    })
  } else {
    out <- sapply(1:size, function(i){
      L0[x$iLtheta] <- theta[i, ]
      cov2cor(tcrossprod(L0))
    })
  }
  dim(out) <- c(p, p, size)
  return(out)
}
