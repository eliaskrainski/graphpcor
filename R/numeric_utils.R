#' Functions for numerical algorithms
#' @name numeric-utils
NULL
#> NULL

#' @describeIn numeric-utils
#' Evaluate the Hessian of the KLD for a `basecor`.
#' @param func model object definition for a correlation matrix.
#' @param x for a `graphpcor` it is the parameter vector,
#' otherwise not used.
#' @param method see [numDeriv::hessian()]
#' @param method.args see [numDeriv::hessian()]
#' @param ... used to pass `ifixed`, an integer
#' vector to indicate model parameters as fixed.
#' If not used, all parameters are treated unknown.
#' @returns matrix with the Hessian
#' @importFrom numDeriv hessian
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
        itheta = func$itheta)
      KLD10(L0 = L0,
            L1 = t(L1))
    },
    x = theta0)
  return(H)
}

#' @describeIn numeric-utils
#' Evaluate the hessian of the KLD for a `basepcor`.
#' @importFrom stats cov2cor
#' @importFrom numDeriv hessian
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
        itheta = func$itheta,
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
#' @importFrom INLAtools is.zero
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
                 itheta = iL, d0 = d0)
  return(hessian(b0, ...))
}
