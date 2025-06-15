#' Build the Cholesky (lower triangular) matrix from theta.
#' @inheritParams pcor
#' @returns matrix with lower triangle as the Cholesky factor
#' of a correlation matrix (if parametrization is 'sap' or 'cpc')
#' or of a precision matrix (if parametrization is 'itp')
#' with diagonal elements as 'd0'.
#' @keywords internal
#' @noRd
theta2L <- function(theta, p, parametrization, itheta, d0) {
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
  if(missing(itheta)) {
    itheta <- which(lower.tri(
      diag(x = rep(1, p), nrow = p, ncol = p)))
  }
  stopifnot(length(theta)==length(itheta))
  stopifnot(all(itheta %in% which(lower.tri(
    diag(x = rep(1, p), nrow = p, ncol = p)))))
  if(missing(d0)) {
    d0 <- d:1
  }
  if(parametrization == "itp") {
    L <- Lprec(
      theta = theta,
      p = p,
      itheta = itheta,
      d0 = d0)
  } else {
    B <- A <- diag(p)
    itheta  <- which(lower.tri(A))
    if(parametrization == 'cpc') {
      A[itheta] <- tanh(theta)
      B[itheta] <- sqrt(1-A[itheta]^2)
    } else {
      theta <- pi/(1+exp(-theta))
      A[itheta] <- cos(theta)
      B[itheta] <- sin(theta)
    }
    for(j in 2:(p-1)) {
      B[, j] <- B[, j] * B[, j-1]
    }
    L <- A * cbind(1, B[, 1:(p-1)])
  }
  attr(L, "parametrization") <- parametrization
  attr(L, "theta") <- theta
  if(parametrization == 'itp') {
    attr(L, "itheta") <- itheta
    attr(L, "d0") <- d0
  }
  return(L)
}
#' Build a correlation matrix from theta
#' @inheritParams pcor
#' @keywords internal
#' @noRd
theta2correl <- function(theta, p, parametrization, itheta, d0) {
  parametrization <- match.arg(
    arg = tolower(parametrization),
    choices = c("cpc", "sap", "itp")
  )
  L <- theta2L(
    theta = theta,
    p = p,
    parametrization = parametrization,
    itheta = itheta,
    d0 = d0
  )
  if(parametrization == 'itp') {
    V <- chol2inv(t(L))
    return(cov2cor(V))
  } else {
    return(tcrossprod(L))
  }
}
#' @describeIn pcor
#' Drawn a random sample correlation matrix.
#' @param lambda numeric to specify the PC-prior parameter,
#' see Simpson et. al. (2017).
#' @param theta0 numeric vector to specify the base model.
#' Drawn \deqn{r \sim \textrm{Exponential}(\lambda)},
#' drawn a vector \deqn{\mathbf{z}} from
#' standard Gaussian distribution.
#' Compute \deqn{\mathbf{s} = \mathbf{z}/\sqrt{\mathbf{z}^{T}\mathbf{z}}}.
#' Compute \deqn{\theta = r\mathbf{H}^{-1/2}\textbf{s}+\theta_0},
#' where \deqn{\mathbf{H}} is the Hessian of the
#' KLD around the base model.
#' @export
rcorrel <- function(lambda, theta0, p, parametrization, itheta, d0) {
  C0 <- theta2correl(theta0, p, parametrization, itheta, d0)
  if(missing(p)) {
    p <- nrow(C0)
  }
  stopifnot(p==nrow(C0))
  m <- p * (p - 1) / 2
  if(missing(lambda))
    lambda <- 0
  if(is.zero(lambda)) {
    theta <- runif(m, 0, pi)
  } else {
    stopifnot(lambda>0.0)
    theta <- INLA::inla.pc.cormat.rtheta(n=1, p, lambda)
  }
  # if(lambda<=sqrt(.Machine$double.eps)) {
  #   g <- rexp(m)
  #   g <- g/sum(g)
  #   theta <- asin(exp(-g))
  #   b <- sample(0:1, size = m, replace = TRUE)
  #   theta <- b * theta + (1-b)*(pi-theta)
  # }
  # if(lambda>sqrt(.Machine$double.eps)) {
  #   theta <- INLA:::inla.pc.cormat.rtheta(n=1, p, lambda)
  # }
  L <- theta2gamma2L(theta, fromR = FALSE)
  R <- tcrossprod(L)
  attr(R, "determinant") <- attr(L, "determinant")
  attr(R, "kld") <- attr(L, "kld")
  return(R)
}
