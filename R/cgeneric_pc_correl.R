#' Build an `cgeneric` object to implement the PC prior,
#' proposed on Simpson et. al. (2007),
#' as an informative prior, see details in [basecor()].
#' @param n integer to define the size of the matrix,
#' same as `p` in [basecor()].
#' @inheritParams basecor
#' @inheritParams cgeneric_graphpcor
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
           itheta,
           iparams,
           ifixed,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           ...) {

    ## collect ... args and check debug
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

    if(missing(base)) {
      if(missing(n)) {
        stop("Please provide 'n' or 'base'!")
      } else {
        stopifnot(n>1)
        warning("Missing base model! Assume zero correlations.")
        base <- diag(n)
      }
    }
    if(inherits(base, "matrix")) {
      if(missing(n)) {
        n <- ncol(base)
      }
    }
    if(is.vector(base)) {
      if(missing(n)) {
        n <- (1 + sqrt(1 + 4 * 2 * length(base))) / 2
      }
    }
    stopifnot(n>1)
    itheta <- p_itheta_fncheck(n, itheta)
    m <- length(itheta)
    if(dotArgs$debug) {
      print(list(n = n, m = m,
                 itheta = itheta))
    }

    iparams <- m_iparams_fncheck(n+m, iparams)
    npars1 <- length(unique(iparams[1:n]))
    npars2 <- length(unique(iparams[n+1:m]))

    basemodel <- basecor(
      base, p = n,
      parametrization = "cpc",
      iparams = iparams[n+1:m]-npars1,
      itheta = itheta
    )
    if(dotArgs$debug) {
      print(basemodel)
    }

    theta0 <- basemodel$theta
    m <- length(theta0)

    if(missing(ifixed)) {
      ifixed <- NULL
    }
    cfixed <- ((1:m) %in% ifixed)

    I0 <- hessian(basemodel, ifixed = ifixed)
    if(dotArgs$debug) {
      cat("Hessian\n")
      print(I0)
    }

    if(is.null(I0)) {
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

    if(dotArgs$debug) {
      print(list(n = n, m = m,
                 iparams = iparams,
                 npars1 = npars1,
                 npars2 = npars2,
                 cfixed = cfixed))
    }

    ## sigma prior
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

    ## update sigma.* with iparams
    sigma.prior.reference <- sigma.prior.reference[iparams[(1:n)]]
    sigma.prior.probability <- sigma.prior.probability[iparams[(1:n)]]
    sigma.fixed <- sigma.fixed[iparams[(1:n)]]
    if(any(is.na(sigma.fixed)))
      stop("Some error in `iparams` or  `sigma.fixed`!")

    if(dotArgs$debug) {
      cat("expanded:\n")
      print(list(sigmaref = sigma.prior.reference,
                 sigmaprob = sigma.prior.probability,
                 sfixed = sigma.fixed))
    }

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
