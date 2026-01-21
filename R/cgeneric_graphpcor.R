#' Build an `cgeneric` for a graph, see [graphpcor()]
#' @description
#' From either a `graph` (see [graph()]) or
#' a square matrix (used as a graph),
#' creates an `cgeneric` (see [INLAtools::cgeneric()])
#' to implement the Penalized Complexity prior using the
#' Kullback-Leibler divergence - KLD from a base graphpcor.
#' @param model  a `graphpcor` (see [graphpcor()]) or
#' a square matrix (to be used as a graph)
#' to define the precision structure of the model.
#' @param lambda the parameter for the exponential prior on
#' the radius of the sphere, see details.
#' @param base numeric vector with length `m`, `m` is the
#' number of edges in the graph, or matrix with the reference
#' correlation model against what the KLD will be evaluated.
#' If it is a vector, a correlation matrix is defined
#' considering the graph model and this vector as
#' the parameters in the lower triangle matrix L.
#' If it is a matrix, it will be checked if the graph model
#' can generates this.
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
#' @param cor.params.fixed numeric vector of length `m`
#' providing the value(s) at which the lower parameter(s)
#' of the L matrix to be fixed and not estimated.
#' NA indicates not fixed and all are set to be estimated by default.
#' Example: with `cor.params.fixed = c(NA, -1, NA, 1)` the first
#' and the third of these parameters will be estimated while
#' the second is fixed and equal to -1 and the forth is fixed
#' and equal to 1. NOTE: `iparams` will be applied here as
#' `cor.params.fixed[iparams[(n+1:m)]-iparams[n]]`, thus the provided
#' examples give `NA -1 -1 NA` and so the second and third low L
#' parameters are fixed to `-1`.
#' @param ... additional arguments that will be passed on to
#' [INLAtools::cgenericBuilder()].
#' @seealso [graphpcor()]
#' @useDynLib graphpcor, .registration = TRUE
#' @returns `cgeneric` object.
cgeneric_graphpcor <-
  function(model,
           lambda,
           base,
           sigma.prior.reference,
           sigma.prior.probability,
           iparams,
           cor.params.fixed,
           ...) {

    ## collect ... args and check debug
    dotArgs <- list(...)
    if(!any(names(dotArgs)=="debug")) {
      dotArgs$debug <- FALSE
    }
    if(inherits(model, "Matrix")) {
      if(dotArgs$debug) {
        cat("Building 'graphpcor' from a 'Matrix'!")
      }
      model <- graphpcor(model)
    }

    ## Matrix structure
    Q0 <- Laplacian(model)
    n <- nrow(Q0)
    stopifnot(n>0)
    if(dotArgs$debug) {
      print(model)
      cat("Laplacian is\n")
      print(Q0)
    }

    if(length(lambda)>1) {
      warning('length(lambda)>1, using lambda[1]!')
    }
    lambda <- as.numeric(lambda[1])
    stopifnot(lambda>0)

    l1 <- t(chol(Q0 + diag(1.0, n, n)))
    qnz <- !is.zero(Q0)
    itheta <- which(qnz & lower.tri(Q0, diag = FALSE))
    qij <- list(
      ii = row(Q0)[itheta],
      jj = col(Q0)[itheta],
      iq = which(Q0!=0))
    qij$ilq <- which(qnz & lower.tri(Q0, diag = TRUE))
    qij$iuq <- which(qnz & upper.tri(Q0, diag = TRUE))
    qij$ilqpac <- which(qnz[lower.tri(Q0, diag = TRUE)])
    ll <- t(chol(Q0 + diag(n)))
    qij$ifil <- setdiff(which(ll!=0), qij$ilq)
    if(dotArgs$debug) {
      print(qij)
    }

    nEdges <- length(qij$ii)
    if(nEdges==0) {
      stop("This graph is trivial, please consider 'iid' model!")
    }
    nnz <- n + nEdges
    nfi <- length(qij$ifil)

    iparams <- m_iparams_fnc(nnz, iparams)
    npars1 <- iparams[n]
    npars2 <- iparams[nnz]-npars1

    if(dotArgs$debug) {
      print(list(n = n, nnz = nnz, nfi = nfi,
                 iparams = iparams,
                 npars1 = npars1, npars2 = npars2))
    }

    if(missing(sigma.prior.reference)) {
       sigma.prior.reference <- rep(1, npars1)
    }
    if(missing(sigma.prior.probability)) {
      sigma.prior.probability <- rep(0, npars1)
    }

    stopifnot(length(sigma.prior.reference) == npars1)
    stopifnot(length(sigma.prior.probability) == npars1)
    stopifnot(all(sigma.prior.reference>0))
    sigma.prior.probability[is.na(sigma.prior.probability)] <- 0
    stopifnot(all(sigma.prior.probability>=0.0))
    stopifnot(all(sigma.prior.probability<=1.0))
    sigma.fixed <- is.zero(sigma.prior.probability) |
      is.zero(1-sigma.prior.probability)
    nUnkSigmas <- npars1 - sum(sigma.fixed)

    ## update sigma.*
    sigma.prior.reference <- sigma.prior.reference[iparams[(1:n)]]
    sigma.prior.probability <- sigma.prior.probability[iparams[(1:n)]]
    sigma.fixed <- sigma.fixed[iparams[(1:n)]]
    if(any(is.na(sigma.fixed)))
      stop("Some error in `iparams` or  `sigma.fixed`!")

    if(dotArgs$debug) {
      print(list(sigmaref = sigma.prior.reference,
                 sigmaprob = sigma.prior.probability,
                 sfixed = sigma.fixed))
    }

    if(missing(cor.params.fixed)) {
      cor.params.fixed <- 1:npars2
    }
    stopifnot(length(cor.params.fixed)==npars2)

    ## expand
    cor.params.fixed <- cor.params.fixed[iparams[n+1:nEdges]-n]
    print(cor.params.fixed)
    if(any(is.na(cor.params.fixed)))
      stop("Some error in `iparams` or  `cor.params.fixed`!")

    ii <- c(1:n, qij$ii)
    jj <- c(1:n, qij$jj)
    ii <- ii[order(jj)]
    jj <- jj[order(jj)]

    iuq <- qij$ilq  ## mem order
    iuqpac <- qij$ilqpac

    ifi <- row(Q0)[qij$ifil]
    jfi <- col(Q0)[qij$ifil]

    if(missing(base)){
      warning("Missing base model! Using 'iid'.")
      base <- rep(0, npars2)
    }

##    I0 <- hessian(model, base, decomposition = "eigen")
    basemodel <- basepcor(
      base,
      p = n,
      itheta = itheta,
      d0 = n:1,
      iparams = iparams[n+1:nEdges]-npars1
    )
    if(dotArgs$debug) {
      cat("base model:\n")
      print(utils::str(basemodel))
    }
    stopifnot(all(dim(basemodel$I0) == c(npars2, npars2)))

    theta0 <- basemodel$theta
    I0 <- basemodel$I0

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
        model = "inla_cgeneric_graphpcor",
        n = as.integer(n),
        debug = as.integer(dotArgs$debug),
        shlib = dotArgs$shlib,
        ne = as.integer(nEdges),
        nfi = as.integer(nfi),
        ii = as.integer(jj-1),
        jj = as.integer(ii-1),
        iuq = as.integer(iuq-1),
        iuqpac = as.integer(iuqpac-1),
        ifi = as.integer(ifi-1),
        jfi = as.integer(jfi-1),
        itheta = as.integer(iparams-1),
        sfixed = as.integer(sigma.fixed),
        lambda = as.numeric(lambda),
        sigmaref = as.numeric(sigma.prior.reference),
        sigmaprob = as.numeric(sigma.prior.probability),
        lconst = as.numeric(attr(I0, "determinant")),
        thetabase = as.numeric(theta0),
        Ihalf = attr(I0, "h.5")
      )
    )
    return(the_model)

  }
