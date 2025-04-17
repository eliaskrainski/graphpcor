#' Build the correlation matrix parametrized from the
#' hypershere decomposition, see details.
#' @rdname correl
#' @param theta numeric vector with length equal n(n-1)/2
#' @param fromR logical indicating if theta is in R.
#' If FALSE, assumes \eqn{\theta[k] \in (0, pi)}.
#' @details
#' The hypershere decomposition, as proposed in
#' Rapisarda, Brigo and Mercurio (2007)
#' consider \eqn{\theta[k] \in [0, \infty], k=1,...,m=n(n-1)/2}
#' compute \eqn{x[k] = pi/(1+exp(-theta[k]))}
#' organize it as a lower triangle of a \eqn{n \times n} matrix
#' \deqn{         | cos(x[i,j])                           ,      j=1}
#' \deqn{B[i,j] = | cos(x[i,j])prod_{k=1}^{j-1}sin(x[i,k]),  2 <= j <= i-1}
#' \deqn{         | prod_{k=1}^{j-1}sin(x[i,k])           ,      j=i}
#' \deqn{         | 0                                     , j+1 <= j <= n }
#' Result
#' \deqn{\gamma[i,j] = -log(sin(x[i,j]))}
#'  \deqn{KLD(R) = \sqrt(2\sum_{i=2}^n\sum_{j=1}^{i-1} \gamma[i,j]}
#' @references
#' Rapisarda, Brigo and Mercurio (2007).
#'   Parameterizing correlations: a geometric interpretation.
#'   IMA Journal of Management Mathematics (2007) 18, 55-73.
#'   <doi 10.1093/imaman/dpl010>
#' @return a correlation matrix
#' @export
theta2correl <- function(theta, fromR = TRUE) {
  tcrossprod(theta2gamma2L(theta, fromR = fromR))
}
#' @describeIn correl
#' Build a lower triangular matrix from a parameter vector.
#' See details.
#' @return Lower triangular n x n matrix
theta2gamma2L <- function(theta, fromR = TRUE) {
  ## see PC-prior paper, section 6.2
  ## return B lower triangular
  m <- length(theta)
  n <- (1 + sqrt(1+8*m))/2
  stopifnot((floor(n)==n) && (ceiling(n)==n) && (n>1))
  if(fromR)
    theta <- pi / (1 + exp(-theta))
  th <- matrix(NA, n, n)
  th[lower.tri(th)] <- theta
  b <- matrix(0, n, n)
  b[1, 1] <- 1
  b[2:n, 1] <- cos(th[2:n, 1])
  su <- sin(th)
  kld <- 0
  for(i in 2:n) {
    kld <- kld + sum(-log(su[i, 1:(i-1)]))
    a <- prod(su[i, 1:(i-1)])
    b[i, i] <- a
    if(i>2) {
      for(j in 2:(i-1)) {
        a <- prod(su[i, 1:(j-1)])
        b[i, j] <- cos(th[i,j]) * a
      }
    }
  }
  attr(b, "determinant") <- prod(diag(b))^2
  attr(b, 'kld') <- sqrt(2 * kld)
  return(b)
}
#' @describeIn correl
#' Drawn a random sample correlation matrix
#' @param p integer to specify the matrix dimension
#' @param lambda numeric as the penalization parameter.
#' If missing it will be assumed equal to zero.
#' The lambda=0 case means no penalization and
#' a random correlation matrix will be drawn.
#' Please see section 6.2 of the PC-prior paper,
#' Simpson et. al. (2017), for details.
#' @references Simspon et. al. (2017).
#' Penalising Model Component Complexity:
#' A Principled, Practical Approach to Constructing Priors.
#'  Statist. Sci. 32(1): 1-28 (February 2017).
#'  <doi: 10.1214/16-STS576>
#' @export
rcorrel <- function(p, lambda) {
  stopifnot(p>1)
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
