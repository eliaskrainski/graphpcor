#' Build an `cgeneric` object to implement the PC prior,
#' proposed on Simpson et. al. (2007),
#' as an informative prior, see details in [basecor()].
#' @param n integer to define the size of the matrix
#' @param lambda numeric (positive), the penalization rate parameter
#' @param base matrix with base correlation matrix,
#' or numeric vector representing the parameters of a base correlation
#' matrix. See [basepcor()] for details.
#' @param sigma.prior.reference numeric vector with length `n`,
#' `n` is the number of nodes (variables) in the graph, as the
#' reference standard deviation to define the PC prior for each
#' marginal variance parameters. If missing, the model will be
#' assumed for a correlation. If a length `n` vector is given
#' and `sigma.prior.probability` is missing, it will be used as
#' known square root of the variances.
#' NOTE: `iparams` will be applied here as
#' `sigma.prior.reference[iparams[1:n]]`.
#' @param sigma.prior.probability numeric vector with length `n`
#' to set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known, as described in
#' `sigma.prior.reference`.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
#' NOTE: `iparams` will be applied here as
#' `sigma.prior.probability[iparams[1:n]]`.
#' @param iparams integer ordered vector with length equals
#' to `n+m` to specify common parameter values. If missing it
#' is assumed `1:(n+m)` and all parameters are assumed distinct.
#' The first `n` indexes the square root of the marginal
#' variances and the remaining indexes the edges parameters.
#' Example: By setting `iparams = c(1,1,2,3, 4,5,5,6)`,
#' the first two standard deviations are common and the
#' second and third edges parameters are common as well,
#' giving 6 unknown parameters in the model.
#' @param ... additional arguments passed on to
#' [INLAtools::cgeneric()].
#' @references
#' Daniel Simpson, H\\aa vard Rue, Andrea Riebler, Thiago G.
#' Martins and Sigrunn H. S\\o rbye (2017).
#' Penalising Model Component Complexity:
#' A Principled, Practical Approach to Constructing Priors.
#' Statistical Science 2017, Vol. 32, No. 1, 1–28.
#' <doi 10.1214/16-STS576>
#' @return a `cgeneric` object, see [cgeneric()] for details.
cgeneric_pc_correl <-
  function(n,
           lambda,
           base,
           sigma.prior.reference,
           sigma.prior.probability,
           iparams,
           ...) {

    if(missing(n) & missing(base))
      stop("Please provide 'n' or 'base'!")
    if(missing(base)) {
      warning("Missing base model! Assume zero correlations.")
      stopifnot(n>1)
      base <- basecor(diag(n))
    } else {
      if(inherits(base, "matrix")) {
        if(missing(n)) {
          n <- ncol(base)
        } else {
          stopifnot(n == ncol(base))
        }
        base <- basecor(base)
      }
      if(is.vector(base)) {
        if(missing(n)) {
          m <- length(base)
          n <- (1 + sqrt(1 + 4 * 2 * m)) / 2
        }
        base <- basecor(base, p = n)
      }
    }
    theta0 <- base$theta
    m <- length(theta0)
    I0 <- hessian(base)

    if(is.null(I0)) {
      I0d <- list(
        logDeterminant = 0,
        sqrt = NULL)
    } else {
      I0d <- dspd(I0)
    }

    stopifnot(lambda>0)

    dotArgs <- list(...)
    if(!any(names(dotArgs)=="debug")) {
      dotArgs$debug <- FALSE
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

    if(length(lambda)>1) {
      warning('length(lambda)>1, using lambda[1]!')
    }
    lambda <- as.numeric(lambda[1])
    stopifnot(lambda>0)

    if(missing(sigma.prior.reference)) {
      sigma.prior.reference <- rep(1, n)
    }
    if(length(sigma.prior.reference)==1) {
      sigma.prior.reference <- rep(sigma.prior.reference, n)
    }
    if(missing(sigma.prior.probability)) {
      sigma.prior.probability <- rep(0, n)
    }
    if(length(sigma.prior.probability)==1) {
      sigma.prior.probability <- rep(sigma.prior.probability, n)
    }
    stopifnot(length(sigma.prior.reference) == n)
    stopifnot(length(sigma.prior.probability) == n)
    stopifnot(all(sigma.prior.reference>0))
    sigma.prior.probability[is.na(sigma.prior.probability)] <- 0
    stopifnot(all(sigma.prior.probability>=0.0))
    stopifnot(all(sigma.prior.probability<=1.0))
    sigma.fixed <- is.zero(sigma.prior.probability) |
      is.zero(1-sigma.prior.probability)

    if(dotArgs$debug) {
      print(list(sigmaref = sigma.prior.reference,
                 sigmaprob = sigma.prior.probability,
                 sfixed = sigma.fixed))
    }

    if(missing(iparams)) {
      iparams <- 1:(n+m)
    } else {
      if(!all(iparams == (1:(n+m))))
        warning("Currently only 1:(n+m)!")
      iparams <- 1:(n+m)
      stopifnot(length(iparams)==(n+m))
      stopifnot(all(iparams %in% (1:(n+m))))
      stopifnot(all(diff(sort(iparams))==1))
    }

    ## update sigmas.prior.*
    sigma.prior.reference <- sigma.prior.reference[iparams[(1:n)]]
    sigma.prior.probability <- sigma.prior.probability[iparams[(1:n)]]
    sigma.fixed <- sigma.fixed[iparams[(1:n)]]
    nSigmas <- iparams[n]
    nLparam <- iparams[n+m]-nSigmas

    the_model <- do.call(
      what = INLAtools::cgenericBuilder,
      args = list(
        model = "inla_cgeneric_pc_correl",
        n = as.integer(n),
        debug = as.logical(dotArgs$debug),
        shlib = dotArgs$shlib,
        itheta = as.integer(iparams -1),
        sfixed = as.integer(sigma.fixed),
        lambda = as.numeric(lambda),
        sigmaref = as.numeric(sigma.prior.reference),
        sigmaprob = as.numeric(sigma.prior.probability),
        lconst = as.numeric(I0d$logDeterminant),
        thetabase = as.numeric(theta0),
        Ihalf = I0d$sqrt
      )
    )

    return(the_model)

}
