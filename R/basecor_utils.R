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

