#' Functions for numerical algorithms
#' @name numeric-utils
NULL
#> NULL

#' @describeIn numeric-utils
#' `hessian`: Evaluate the Hessian of the KLD for a `basecor`.
#' `spd`: decompose a positive definite matrix,
#' and compute useful elements out of that.
#' @param M square matrix.
#' @param decomposition character to inform
#' which decomposition is to be applied to the
#' hessian. The options are "eigen", "svd" and "chol".
#' Default is "svd".
#' @returns `dspd` returns a list with the decomposition elements,
#'  "logDeterminant" (of the original matrix),
#' "sqrt" (its 'square root') and
#' "sqrtInv" (its inverse 'square root').
#' `hessian` returns the Hessian
dspd <- function(M, decomposition = "svd") {
  decomposition <- match.arg(
    tolower(decomposition),
    c("svd", "eigen", "chol"))
  stopifnot(nrow(M)==ncol(M))
  p <- ncol(M)
  ## next bit follows mvtnorm:::rmvnorm()
  t0 <- sqrt(.Machine$double.eps)
  if(decomposition == "eigen") {
    xd <- eigen(M)
    tol1 <- t0 * abs(xd$values[1])
    if(!all(xd$values >= tol1))
      warning("'M' is numerically not positive semidefinite")
    xd$logDeterminant <- sum(log(xd$values))
    s <- sqrt(pmax(xd$values, 0.0))
    xd$sqrt <- t(xd$vectors %*% (t(xd$vectors) * s))
    xd$sqrtInv <- t(xd$vectors %*% (t(xd$vectors) / s))
  }
  if(decomposition == "svd") {
    xd <- svd(M)
    tol1 <- t0 * abs(xd$d[1])
    if(any(xd$d<tol1))
      warning("'M' is numerically not positive semidefinite")
    xd$logDeterminant <- sum(log(xd$d))
    s <- sqrt(pmax(xd$d, 0.0))
    xd$sqrt <- t(xd$v %*% (t(xd$u) * s))
    xd$sqrtInv <- t(xd$v %*% (t(xd$u) / s))
  }
  if(decomposition == "chol") {
    xd <- chol(M, pivot = TRUE)
    if(any(diag(xd)<t0))
      warning("'M' is numerically not positive semidefinite")
    tol1 <- pmin(t0 * 10, diag(xd))
    xd$logDeterminant <- sum(diag(xd))*2.0
    x$sqrt <- matrix(xd[, order(attr(xd, "pivot")), ], p)
    xn <- chol2inv(chol(M))
    xn.5 <- chol(xn, pivot = TRUE)
    xd$sqrtInv <- matrix(xn.5[, order(attr(xn.5, "pivot"))], p)
  }
  return(xd)
}
#' @describeIn numeric-utils
#' Evaluate the Hessian for a `basecor`.
#' @param func model object definition for a correlation matrix.
#' @param x for a `graphpcor` it is the parameter vector,
#' otherwise not used.
#' @param method see [numDeriv::hessian()]
#' @param method.args see [numDeriv::hessian()]
#' @param ... used to pass `ifixed`, an integer
#' vector to indicate model parameters as fixed.
#' If not used, all parameters are treated unknown.
#' @returns matrix with the Hessian
#' @export
hessian.basecor <- function(
    func,
    x,
    method = "Richardson",
    method.args = list(),
    ...) {

  ifixed <- list(...)$ifixed
  m <- length(func$theta)
  nunk <- m - length(ifixed)
  if(nunk==0) return(NULL)

  iunknown <- setdiff(1:m, ifixed)
  theta0 <- func$theta[iunknown]

  ## the KLD10() uses upper triangle Cholesky
  L0 <- t(func$L)
  H <- hessian(
    func = function(x) {
      th <- func$theta
      th[iunknown] <- x
      L1 <- cholcor(
        theta = th[func$iparams],
        p = func$p,
        parametrization = func$parametrization,
        iLtheta = func$iLtheta)
      KLD10(L0 = L0,
            L1 = t(L1))
    },
    x = theta0)
  return(H)
}

#' @describeIn numeric-utils
#' Evaluate the hessian of the KLD for a `basepcor`.
#' @export
hessian.basepcor <- function(
    func,
    x,
    method = "Richardson",
    method.args = list(),
    ...) {

  ifixed <- list(...)$ifixed
  m <- length(func$theta)
  nunk <- m - length(ifixed)
  if(nunk==0) return(NULL)

  iunknown <- setdiff(1:m, ifixed)
  theta0 <- func$theta[iunknown]

  L0 <- func$L ## the correlation's (lower) Cholesky
  H <- hessian(
    func = function(x) {
      th <- func$theta
      th[iunknown] <- x
      L1Q0 <- Lprec0(
        theta = th[func$iparams],
        p = func$p,
        iLtheta = func$iLtheta,
        d0 = func$d0)
      C1 <- cov2cor(chol2inv(t(L1Q0)))
      return(KLD10(L0 = L0, C1 = C1))
    },
    x = theta0)
  return(H)
}

#' @describeIn numeric-utils
#' Evaluate the hessian of the KLD for a `graphpcor`
#' correlation model around a base model.
#' @export
hessian.graphpcor <- function(
    func,
    x,
    method = "Richardson",
    method.args = list(),
    ...) {

  d <- dim(func)
  iL <- which(lower.tri(diag(d[1])) &
                !is.zero(Laplacian(func)))
  dotArgs <- list(...)
  if(any(names(dotArgs)=="d0")) {
    d0 <- dotArgs$d0
  } else {
    d0 <- d[1]:1
  }
  b0 <- basepcor(base = x, p = d[1],
                 iLtheta = iL, d0 = d0)
  return(hessian(b0, ...))
}
