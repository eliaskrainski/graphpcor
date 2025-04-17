#' The LKJ density for a correlation matrix
#' @param R correlation matrix
#' @param eta numeric, the prior parameter
#' @param log logical indicating if the log of the density
#' is to be returned, default = FALSE
#' @return numeric as the (log) density
#' @export
dLKJ <- function(R, eta, log = FALSE) {
  lR <- chol(R)
  ldR <- 2*sum(log(diag(lR)))
  d <- ncol(R)
  k <- 1:(d-1)
  lbk <- lbeta(eta + (d-k-1)/2,
               eta + (d-k-1)/2) * (d-k)
  p2 <- sum((2*eta -2 + d - k)*(d-k))
  o <- sum(lbk) + p2*log(2) + (eta-1)*ldR
  if(!log)
    o <- exp(o)
  return(o)
}
#' Build an `inla.cgeneric` object to implement the
#' LKG prior for the correlation matrix.
#' @param n integer to define the size of the matrix
#' @param eta numeric greater than 1, the parameter
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param libpath string, default is NULL, with the path to the shared object.
#' @details
#' The parametrization uses the
#' hypershere decomposition, as proposed in
#' Rapisarda, Brigo and Mercurio (2007).
#' consider \eqn{\theta[k] \in [0, \infty], k=1,...,m=n(n-1)/2}
#' from \eqn{\theta[k] \in [0, \infty], k=1,...,m=n(n-1)/2}
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
#' @return a `inla.cgeneric`, [cgeneric()] object.
cgeneric_LKJ <-
  function(n,
           eta,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    if(is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("graphpcor")
      } else {
        libpath <- system.file("libs", package = "graphpcor")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "graphpcor.dll")
        } else {
          libpath <- file.path(libpath, "graphpcor.so")
        }
      }
    }

    stopifnot(n>1)
    stopifnot(eta>0)

    k <- 1:(n-1)
    lc <- sum((2*eta-2+n-k)*(n-k))*log(2) +
      sum(lbeta(eta + (n-k-1)/2,
                eta + (n-k-1)/2)*(n-k))

    if(debug) {
      cat('log C', lc, '\n')
    }

    cmodel = "inla_cgeneric_LKJ"

    the_model <- list(
      f = list(
        model = "cgeneric",
        n = as.integer(n),
        cgeneric = list(
          model = cmodel,
          shlib = libpath,
          n = as.integer(n),
          debug = as.logical(debug),
          data = list(
            ints = list(
              n = as.integer(n),
              debug = as.integer(debug)
            ),
            doubles = list(
              eta = as.double(eta),
              lc = as.double(lc)
            ),
            characters = list(
              model = cmodel,
              shlib = libpath
            ),
            matrices = list(
              ),
            smatrices = list(
              )
            )
          )
        )
      )

    class(the_model) <- "inla.cgeneric"
    class(the_model$f$cgeneric) <- "inla.cgeneric"

    return(the_model)

}
