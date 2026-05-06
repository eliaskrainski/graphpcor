#' Build an `cgeneric` for a correlation matrix PC-prior.
#' @param n integer to define the size of the matrix,
#' same as `p` in [basecor()].
#' @param base numeric vector, matrix  or `basecor` to define the base
#' correlation model. See [basecor()] for details.
#' If the output of a [basecor()] is provided,
#' `iLtheta` and `iparams` (for the correlation parameters)
#' will be considered from this.
#' @param iLtheta integer vector to specify the (vectorized) position
#' where 'theta' will be placed in the (lower triangle)  Cholesky
#' factorization of the correlation matrix.
#' @details
#' The parametrization is set as in [basecor()] and the base
#' is used to define an informative prior, as derived in
#' the pcmultivariate vignette.
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
           base,
           iLtheta,
           iparams,
           ...) {

    if((!missing(base)) && inherits(base, "basecor")) {
      return(cgeneric(
        model = base, ...))
    }

    dotArgs <- list(...)
    if(!any(names(dotArgs)=="debug")) {
      dotArgs$debug <- FALSE
    }

    ## check 'n' and 'iLtheta'
    if((!missing(base)) && inherits(base, "matrix"))
      n <- ncol(base)
    iLtheta <- p_iLtheta_fncheck(n, iLtheta)
    n <- attr(iLtheta, "p")
    m <- length(iLtheta)

    ## check 'iparams'
    iparams <- m_iparams_fncheck(n + m, iparams)
    npars1 <- length(unique(iparams[1:n]))
    stopifnot(npars1>0)
    npars2 <- length(unique(iparams[-(1:n)]))
    stopifnot(npars2>0)

    if(dotArgs$debug) {
      print(list(n = n, iparams = iparams,
                 npars1 = npars1, npars2 = npars2))
    }

    if(missing(base) || is.null(base)) {
      warning("Missing base model! Assume zero correlations.")
      base <- rep(0.0, npars2)
    } else {
      if(is.vector(base)) {
        stopifnot(length(base)==npars2)
      }
    }

    if(!inherits(base, "basecor")) {
      base <- basecor(
        base = base,
        p = n,
        iparams = iparams[-(1:n)]-iparams[n],
        iLtheta = iLtheta
      )
    }

    return(cgeneric(
      model = base,
      iparams = iparams,
      ...))
}
#' @describeIn cgeneric_pc_correl
#' Build a `cgeneric` for a `basecor`.
#' @param model a `basecor` object.
#' @inheritParams cgeneric_graphpcor
#' @export
cgeneric.basecor <-
  function(model,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           iparams,
           cfixed,
           ...) {

  stopifnot(inherits(model, 'basecor'))
  n <- model$p
  stopifnot(n>1)

  stopifnot(length(model$iLtheta)>0)
  m1 <- length(model$iparams)
  stopifnot(length(model$iLtheta)==m1)

  ## collect ... args and check debug
  dotArgs <- list(...)
  if(!any(names(dotArgs)=="debug")) {
    dotArgs$debug <- FALSE
  }

  if(is.null(dotArgs$useINLAprecomp) ||
     dotArgs$useINLAprecomp) {
    vs <- "26.03.19"
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

  if(length(lambda)>1) {
    warning('length(lambda)>1, using lambda[1]!')
  }
  lambda <- as.numeric(lambda[1])
  stopifnot(lambda>0)

  ## note: iparams for sigma and correl params
  m <- length(model$iparams)
  iparams <- m_iparams_fncheck(n+m, iparams)
  stopifnot(all(model$iparams == (iparams[-(1:n)]-iparams[n])))
  npars1 <- length(unique(iparams[1:n]))
  stopifnot(npars1>0)
  npars2 <- length(unique(iparams[-(1:n)]))
  stopifnot(npars2>0)
  if(dotArgs$debug) {
    print(list(n = n, iparams = iparams,
               npars1 = npars1, npars2 = npars2))
  }

  if(dotArgs$debug) {
    print(model)
  }
  theta0 <- model$theta

  if(missing(cfixed)) {
    cfixed <- NULL
  }
  cfixed <- ((1:npars2) %in% cfixed)
  I0 <- hessian(model, ifixed = which(cfixed))
  if(dotArgs$debug) {
    cat("Hessian\n")
    print(I0)
  }

  if(is.null(I0)) {
    stop("This is not suppose to happen!")
    I0d <- list(
      logDeterminant = 0,
      sqrt = NULL)
  } else {
    I0d <- dspd(I0)
  }
  if(dotArgs$debug) {
    cat("Hessian decomposition\n")
    print(I0d)
  }

  ## sigma prior
  if(is.na(packageCheck(
    name = "INLAtools",
    minimum_version = "0.1.1.904"
  ))) {
    stop("Please update INLAtools")
  }
  pcSigmas <- INLAtools::pcParamCheck(
    npars = npars1,
    reference = sigma.prior.reference,
    probability = sigma.prior.probability
  )
  if(dotArgs$debug) {
    print(pcSigmas)
  }

  iLvec <- pmatch(model$iLtheta,
                  which(lower.tri(matrix(0, n, n))))

  the_model <- do.call(
    what = INLAtools::cgenericBuilder,
    args = list(
      model = "inla_cgeneric_pc_correl",
      n = as.integer(n),                              ## i0
      debug = as.logical(dotArgs$debug),              ## i1
      shlib = dotArgs$shlib,
      iLtheta = as.integer(iLvec-1),                  ## i2
      ifixed = as.integer(c(pcSigmas$fixed, cfixed)), ## i3
      iparams = as.integer(iparams-1),                ## i4
      lambda = as.numeric(lambda),                   ## d0
      sigmaref = as.numeric(pcSigmas$reference),     ## d1
      sigmaprob = as.numeric(pcSigmas$probability),  ## d2
      lconst = as.numeric(I0d$logDeterminant * 0.5), ## d3
      thetabase = as.numeric(theta0),                ## d4
      Ihalf = I0d$sqrt
    )
  )

  return(the_model)

}
