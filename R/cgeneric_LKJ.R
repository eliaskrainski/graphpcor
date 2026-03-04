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
  o <- (eta-1)*ldR - sum(lbk) - p2*log(2)
  if(!log)
    o <- exp(o)
  return(o)
}
#' Build an `cgeneric` model for
#' the LKG prior on correlation matrix.
#' @param n integer to define the size of the matrix
#' @param eta numeric greater than 1, the parameter
#' @param sigma.prior.reference numeric vector with length `n`,
#' `n` is the number of nodes (variables) in the graph, as the
#' reference standard deviation to define the PC prior for each
#' marginal variance parameters. If missing, the model will be
#' assumed for a correlation. If a length `n` vector is given
#' and `sigma.prior.reference` is missing, it will be used as
#' known square root of the variances.
#' @param sigma.prior.probability numeric vector with length `n`
#' to set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known, as described in
#' `sigma.prior.reference`.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
#' @param ... additional arguments passed to [INLAtools::cgeneric()].
#' @seealso It uses the Cannonical Partial Correlation (CPC),
#' parametrization, see [basecor()] for details.
#' @return a `cgeneric` object, see [INLAtools::cgeneric()] for details.
cgeneric_LKJ <-
  function(n,
           eta,
           sigma.prior.reference = rep(1, n),
           sigma.prior.probability = rep(NA, n),
           ...) {

    stopifnot(n>1)
    stopifnot(eta>0)

    dotArgs <- list(...)
    if(!any(names(dotArgs)=="debug")) {
      dotArgs$debug <- FALSE
    }

    scheck <- pcSigmasCheck(
      nsigmas = n,
      sigma.prior.reference = sigma.prior.reference,
      sigma.prior.probability = sigma.prior.probability
    )
    if(dotArgs$debug) {
      print(scheck)
    }

    k <- 1:(n-1)
    lc <- sum((2*eta-2+n-k)*(n-k))*log(2) +
      sum(lbeta(eta + (n-k-1)/2,
                eta + (n-k-1)/2)*(n-k))

    if(dotArgs$debug) {
      cat("log C = ", lc, "\n")
    }

    if(is.null(dotArgs$shlib)) {
      if(dotArgs$debug){
        cat("searching shlib...\n")
      }
      dotArgs$shlib <- do.call(
        what = INLAtools::cgeneric_shlib,
        args = c(list(package = "graphpcor"),
                    dotArgs))
    }

    the_model <- do.call(
      what = INLAtools::cgenericBuilder,
      args = list(
        model = "inla_cgeneric_LKJ",
        n = as.integer(n),
        debug = dotArgs$debug,
        eta = as.double(eta),
        lc = as.double(lc),
        shlib = dotArgs$shlib,
        sfixed = as.integer(scheck$sigma.fixed),
        sigmaref = as.numeric(scheck$sigma.prior.reference),
        sigmaprob = as.numeric(scheck$sigma.prior.probability)
      )
    )

    return(the_model)

}
