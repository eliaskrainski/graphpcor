#' Internal functions.
#' @name internal-utils
NULL
#> NULL

#' @describeIn internal-utils
#' Function to deal with `p` and `iLtheta`
#' @inheritParams basepcor
p_iLtheta_fncheck <- function(p, iLtheta) {
  if(missing(iLtheta)) {
    if(missing(p))
      stop("Please provide 'p' or 'iLtheta'!")
    p <- as.integer(p[1])
    iLtheta <- which(lower.tri(matrix(0, p, p)))
  } else {
    if(inherits(iLtheta, "graphpcor")) {
      Q1 <- as.matrix(Laplacian(iLtheta))
      p <- ncol(Q1)
      iLtheta <- which(lower.tri(Q1) & (Q1 != 0.0))
    } else {
      if(missing(p)) {
        stop("Please provide 'p'!")
      }
    }
    if(is.null(iLtheta)) {
      iLtheta <- which(lower.tri(matrix(0, p, p)))
    }
  }
  il0 <- which(lower.tri(matrix(0, p, p)))
  stopifnot(all(iLtheta %in% il0))
  attr(iLtheta, 'p') <- as.integer(p)
  return(iLtheta)
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
#' Check the PC-prior arguments for `sigma`.
#' @param nsigmas number of parameters.
#' @param sigma.prior.reference numeric vector to set the reference
#' for each standard deviation parameter for its PC-prior.
#' @param sigma.prior.probability numeric vector with to
#' set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
pcSigmasCheck <- function(nsigmas,
                         sigma.prior.reference,
                         sigma.prior.probability) {
  if(missing(sigma.prior.reference)) {
    sigma.prior.reference <- rep(1, nsigmas)
  }
  if(length(sigma.prior.reference)==1) {
    sigma.prior.reference <- rep(sigma.prior.reference, nsigmas)
  }
  if(missing(sigma.prior.probability)) {
    sigma.prior.probability <- rep(0, nsigmas)
  }
  if(length(sigma.prior.probability)==1) {
    sigma.prior.probability <- rep(sigma.prior.probability, nsigmas)
  }

  stopifnot(length(sigma.prior.reference) == nsigmas)
  stopifnot(length(sigma.prior.probability) == nsigmas)
  stopifnot(all(sigma.prior.reference>0))
  sigma.prior.probability[is.na(sigma.prior.probability)] <- 0
  stopifnot(all(sigma.prior.probability>=0.0))
  stopifnot(all(sigma.prior.probability<=1.0))
  sigma.fixed <- is.zero(sigma.prior.probability) |
    is.zero(1-sigma.prior.probability)
  return(list(
    sigma.prior.reference = sigma.prior.reference,
    sigma.prior.probability = sigma.prior.probability,
    sigma.fixed = sigma.fixed
  ))
}
