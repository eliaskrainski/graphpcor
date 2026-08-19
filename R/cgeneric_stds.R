#' Build an `cgeneric` for an iid model with different variances.
#' @description
#' The precision/covariance is diagonal with each entry a
#' different parameter.
#' @param n the model dimension.
#' @param sigma.prior.reference numeric vector to set the reference
#' for each standard deviation parameter for its PC-prior.
#' If missing, the model will be assumed as for a correlation.
#' Note: `iparams` will be applied here as
#' `sigma.prior.reference[iparams[1:n]]`.
#' @param sigma.prior.probability numeric vector with to
#' set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
#' Note: `iparams` will be applied here as
#' `sigma.prior.probability[iparams[1:n]]`.
#' @param iparams integer vector with length equal `n+m`,
#' where `m` is the number of correlation parameters,
#' to identify (possible) common parameters in the model.
#' Default is `1:(n+m)`.
#' Note: `c(1,2,1)` is allowed, but `c(2,1,2)` is not.
#' @param ... additional arguments passed on to
#' [INLAtools::cgeneric()], such as `debug`,
#' `shlib` and `useINLAprecomp`.
#' @details
#' The parametrization is set as theta_i=log(sigma_i).
#' @returns `cgeneric` object.
cgeneric_stds <-
  function(n,
           sigma.prior.reference,
           sigma.prior.probability,
           iparams,
           ...) {

    stopifnot(n>0)

    ## collect ... args and check debug
    dotArgs <- list(...)
    if(!any(names(dotArgs)=="debug")) {
      dotArgs$debug <- FALSE
    }

    if(is.null(dotArgs$useINLAprecomp) ||
       dotArgs$useINLAprecomp) {
      vs <- "26.06.10"
      ivs <- packageCheck(
        name = "INLA",
        minimum_version = vs) >= vs
      if(is.na(ivs)) {
        warning("Update INLA to use 'useINLAprecomp = TRUE'\n")
        dotArgs$useINLAprecomp = FALSE
      }
    }

    if(is.null(dotArgs$shlib)) {
      if(dotArgs$debug){
        cat("searching shlib...\n")
      }
      dotArgs$shlib <-
        INLAtools::cgeneric_shlib_path(
          package = "graphpcor",
          useINLAprecomp = dotArgs$useINLAprecomp,
          debug = dotArgs$debug
          )
    }
    ## sigma prior
    if(is.na(packageCheck(
      name = "INLAtools",
      minimum_version = "0.1.1.904"
    ))) {
      stop("Please update INLAtools")
    }

    ## check iparams
    if(missing(iparams) || is.null(iparams)) {
      iparams <- 1:n
    }
    iparams <- m_iparams_fncheck(iparams=iparams)
    npars <- attr(iparams, "m")

    pcSigmas <- INLAtools::pcParamCheck(
      npars = npars,
      reference = sigma.prior.reference,
      probability = sigma.prior.probability
    )
    if(dotArgs$debug) {
      print(pcSigmas)
    }

    the_model <- do.call(
      what = INLAtools::cgenericBuilder,
      args = list(
        model = "inla_cgeneric_stds",
        n = as.integer(n),                 ## i0
        debug = as.integer(dotArgs$debug), ## i1
        shlib = dotArgs$shlib,
        iparams = as.integer(iparams-1),          ## i2
        ifixed = as.integer(pcSigmas$fixed),       ## i3
        sigmaref = as.numeric(pcSigmas$reference),    ## d0
        sigmaprob = as.numeric(pcSigmas$probability)  ## d1
      )
    )
    return(the_model)

  }
