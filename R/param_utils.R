#' Internal function to interpret theta for precision
#' @rdname param-utils
#' @param theta numeric vector of length m = p(p-1)/2
#' @param p integer with the resulting matrix dimension
interpret.ptheta <- function(theta, p) {
  ii <- 1:p
  ij <- p + 1:(p * (p-1)/2)
  L.prec <- diag(exp(theta[ii]))
  L.prec[lower.tri(L.prec)] <- theta[ij]
  ret <- list(mprec = tcrossprod(L.prec))
  ret$covar <- chol2inv(t(L.prec))
  ret$sd <- sqrt(diag(ret$covar))
  ret$correl <- t(ret$covar/ret$sd)/ret$sd
  return(ret)
}
#' @describeIn param-utils
#' Internal function to interpret theta for covariance
interpret.vtheta <- function(theta, p) {
  ii <- 1:p
  ij <- p + 1:(p * (p-1)/2)
  Lcov <- diag(exp(-theta[ii]))
  Lcov[lower.tri(Lcov)] <- theta[ij]
  ret <- list(mprec = chol2inv(t(Lcov)))
  ret$covar <- tcrossprod(Lcov)
  ret$sd <- sqrt(diag(ret$covar))
  ret$correl <- t(ret$covar/ret$sd)/ret$sd
  return(ret)
}
#' @describeIn param-utils
#' From theta to the lower triangle such that Q = LL'.
#' @param theta numeric vector of length m with
#' the parameters for the lower part of L.
#' If Q is dense, then `m = p(p-1)/2`, else
#' `m = length(ilowerL)`.
#' @param p numeric giving the dimention of Q.
#' If missing, `p = (1+sqrt(1+8*length(theta)))`
#' and Q is assumed to be dense.
#' @param ilowerL numeric vector with index for the
#' part of L to be filled with `theta`.
#' Default is missing and Q is assumed to be dense.
#' @examples
#' theta1 <- c(1, -1, -2)
#' graphpcor:::theta2L(theta1)
#' graphpcor:::theta2L(theta1, 4, c(2,4,7,12))
theta2L <- function(theta, p, ilowerL) {
  ### imput theta output L such Q = LL'
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
  }
  return(L)
}
#' @describeIn param-utils
#' From theta to a correlation matrix
#' using the 'theta2L' parametrization
#' @importFrom stats cov2cor
#' @examples
#' graphpcor:::theta2C(c(0,0,0))
#' graphpcor:::theta2C(c(1,1,1))
theta2C <- function(theta) {
### imput theta output C
    V <- chol2inv(t(theta2L(theta)))
    return(cov2cor(V))
}
#' @describeIn param-utils
#' Compute the KLD with respect to a base model
#' @param C1 is a correlation matrix.
#' @param C0 is a correlation matrix of the base model.
#' @details
#' compute C1 using 'theta2C' on theta  with
#'  \deqn{KLD = 0.5( tr(C0^{-1}C1) -p + ... - log(|C1|) + log(|C0|) )}
#' @examples
#' cc0 <- graphpcor:::theta2C(c(-3,-3,-3))
#' cc1 <- graphpcor:::theta2C(c(-1,-1,-1))
#' graphpcor:::KLD10(cc1, cc0)
KLD10 <- function(C1, C0) {
### imput C1, C0 ouptut KLD
    p <- nrow(C1)
    l1 <- chol(C1)
    hld1 <- sum(log(diag(l1)))
    if(missing(C0)) {
        warning("Missing base model!")
        C0 <- diag(rep(1, p), p, p)
    }
    l0 <- chol(C0)
    hld0 <- sum(log(diag(l0)))
    tr <- sum(diag(chol2inv(l0) %*% C1))
    0.5*(tr -p) + hld0 - hld1
}
#' @describeIn param-utils
#' Compute the hessian, its svd and some elements
#' @examples
#' graphpcor:::theta2H(c(-1,-1,-1))
theta2H <- function(theta) {
### imput theta output H, H^0.5, svd(H)
    if(missing(theta)) stop('Please provide "theta"!')
    C0 <- theta2C(theta)
    H <- hessian(function(x) KLD10(theta2C(x), C0),
                 theta)
    sv <- svd(H)
    if(any(sv$d<sqrt(.Machine$double.eps) * abs(sv$d[1])))
        warning("'H' is numerically not positive definite")
    h.5 <- t(sv$v %*% (t(sv$u) * sqrt(sv$d)))
    hneg.5 <- t(sv$v %*% (t(sv$u) / sqrt(sv$d)))
    stopifnot(all.equal(H, tcrossprod(h.5)))
    list(H = H,
         h.5 = h.5,
         hneg.5 = hneg.5,
         svdH = sv)
}
#' @describeIn param-utils
#' Transform from spherical coordinates to Euclidian
#' @param rphi numeric vector where the first element
#' is the radius and the remaining are the angles
#' @details
#' see [N-sphere/Euclidian](https://en.wikipedia.org/wiki/N-sphere)
rphi2x <- function(rphi) {
### to convert from \{r, \phi_1, ..., \phi_{m-1} \} into x_i \in \Re
### see https://en.wikipedia.org/wiki/N-sphere
    co <- cos(rphi[-1])
    si <- cumprod(sin(rphi[-1]))
    x <- rphi[1] * c(co, 1) * c(1, si)
    return(x)
}
#' @describeIn param-utils
#' Tranform from Euclidian coordinates to spherical
#' @param x parameters in the Euclidian space to be converted
x2rphi <- function(x) {
### to convert from x_i \in \Re into \{r, \phi_1, ..., \phi_{m-1} \}
### see https://en.wikipedia.org/wiki/N-sphere
### NOTE: from x to phi it may give phi[m-1]<0, if so add 2*pi
    m <- length(x)
    if(m>1) {
      phi <- numeric(m-1)
      phi[m-1] <- atan2(x[m], x[m-1])
      if(phi[m-1]<0)
        phi[m-1] <- phi[m-1] + 2*pi
      r2 <- x[m]^2 + x[m-1]^2
      if(m>2) {
        for(i in (m-2):1) {
          phi[i] <- atan2(sqrt(r2), x[i])
          r2 <- r2 + x[i]^2
        }
      }
    } else {
      r2 <- x^2
      phi <- NULL
    }
    return(c(sqrt(r2), phi))
}
#' @describeIn param-utils
#' Drawn samples from the PC-prior for correlation
#' @param n integer to define the size of the correlation matrix
#' @param lambda numeric as the parameter for the
#' Exponential distribution of the radius
#' @param R scaling matrix (square root of the Hessian
#' around the base model)
#' @param theta.base numeric vector of the base model
#' @importFrom stats runif
#' @importFrom stats rexp
rtheta <- function(n, lambda=1, R, theta.base) {
    m <- n*(n-1)/2
    r <- rexp(1, lambda) ## radial coordinate
    if(m>1) {
      phi <- numeric(m-1)
      if(m>2)
        phi[1:(m-2)] <- runif(m-2, 0, pi)  ## m-2 angles
      phi[m-1] <- runif(1, 0, 2*pi)  ## last angle

    } else {
      phi <- NULL
    }
    drop((rphi2x(c(r, phi)) %*% R) + theta.base)
}
#' @describeIn param-utils
#' PC-prior density for the correlation matrix
#' @param H.elements list output of theta2H
dtheta <- function(theta, lambda, theta.base, H.elements) {
    ## log(\lambda)-r\lamda -log(2)-(m-1)log(\pi) +log(|J1|) +log(|J2|)
    ## log(|J1|) = {m-1}*log(r) + \sum_{i=1:m-1} (m-i-1)*log(sin(\phi_i))
    ## log(|J2|) is numerically evaluated as log(abs(det(H)))
    stopifnot(length(theta)==length(theta.base))
    m <- length(theta)
    ld2 <- 2*sum(log(abs(H.elements$svd$d))) ## log(abs(det(H)))
    x <- drop(H.elements$hneg.5 %*%theta) - theta.base
    rphi <- x2rphi(x)
    stopifnot(all(rphi>0))
    ld1 <- (m-1) * log(rphi[1])
    if(m>2) {
      ii <- 1:(m-2) ### do not include the last one!
      ld1 <- ld1 + sum(log(sin(rphi[ii+1])) * (m-ii-1))
    }
    ld <- log(lambda) -rphi[1]*lambda -log(2) - (m-1)*log(pi)
    return(ld + ld1 + ld2)
}
