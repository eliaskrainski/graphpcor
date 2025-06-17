#' Build an `inla.cgeneric` for a graph, see [graphpcor()]
#' @description
#' From either a `graph` (see [graph()]) or
#' a square matrix (used as a graph),
#' creates an `inla.cgeneric` (see [cgeneric()])
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
#' and `sigma.prior.reference` is missing, it will be used as
#' known square root of the variances.
#' NOTE: `params.id` will be applied here as
#' `sigma.prior.reference[params.id[1:n]]`.
#' @param sigma.prior.probability numeric vector with length `n`
#' to set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known, as described in
#' `sigma.prior.reference`.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
#' NOTE: `params.id` will be applied here as
#' `sigma.prior.probability[params.id[1:n]]`.
#' @param params.id integer ordered vector with length equals
#' to `n+m` to specify common parameter values. If missing it
#' is assumed `1:(n+m)` and all parameters are assumed distinct.
#' The first `n` indexes the square root of the marginal
#' variances and the remaining indexes the edges parameters.
#' Example: By setting `params.id = c(1,1,2,3, 4,5,5,6)`,
#' the first two standard deviations are common and the
#' second and third edges parameters are common as well,
#' giving 6 unknown parameters in the model.
#' @param low.params.fixed numeric vector of length `m`
#' providing the value(s) at which the lower parameter(s)
#' of the L matrix to be fixed and not estimated.
#' NA indicates not fixed and all are set to be estimated by default.
#' Example: with `low.params.fixed = c(NA, -1, NA, 1)` the first
#' and the third of these parameters will be estimated while
#' the second is fixed and equal to -1 and the forth is fixed
#' and equal to 1. NOTE: `params.id` will be applied here as
#' `low.params.fixed[params.id[(n+1:m)]-n+1]`, thus the provided
#' examples give `NA -1 -1 NA` and so the second and third low L
#' parameters are fixed to `-1`.
#' @param ... additional arguments that will be passed on to
#' [INLAtools::cgeneric()] such as
#' `debug`, `useINLAprecomp` and `libpath`.
#' @returns `cgeneric`/`inla.cgeneric` object.
#' @useDynLib graphpcor, .registration = TRUE
cgeneric_graphpcor <-
  function(model,
           lambda,
           base,
           sigma.prior.reference,
           sigma.prior.probability,
           params.id,
           low.params.fixed,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    if (is.null(libpath)) {
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

    if(inherits(model, "matrix")) {
      model <- graphpcor(model)
    }
    Q0 <- Laplacian(model)
    n <- nrow(Q0)
    stopifnot(n>0)
    if(debug>99) {
      print(model)
      cat("Laplacian is\n")
      print(Q0)
    }
    stopifnot(all(lambda>0))
    if(length(lambda)>1) {
      warning('length(lambda)>1, using lambda[1]!')
      lambda <- as.numeric(lambda[1])
    }
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

    if(debug) {
      print(list(sigmaref = sigma.prior.reference,
                 sigmaprob = sigma.prior.probability,
                 sfixed = sigma.fixed))
    }

    l1 <- t(chol(Q0 + diag(1.0, n, n)))
    qnz <- !is.zero(Q0)
    qij <- list(
      ii = row(Q0)[qnz & lower.tri(Q0, diag = FALSE)],
      jj = col(Q0)[qnz & lower.tri(Q0, diag = FALSE)],
      iq = which(Q0!=0))
    qij$ilq <- which(qnz & lower.tri(Q0, diag = TRUE))
    qij$iuq <- which(qnz & upper.tri(Q0, diag = TRUE))
    qij$ilqpac <- which(qnz[lower.tri(Q0, diag = TRUE)])
    ll <- t(chol(Q0 + diag(n)))
    qij$ifil <- setdiff(which(ll!=0), qij$ilq)
    if(debug>99) {
      print(qij)
    }

    nEdges <- length(qij$ii)
    nnz <- n + nEdges
    nfi <- length(qij$ifil)

    if(missing(params.id)) {
      params.id <- 1:nnz
    } else {
      stopifnot(length(params.id)==nnz)
      stopifnot(all(params.id %in% (1:nnz)))
      stopifnot(all(diff(sort(params.id))>0))
    }
    ## update sigmas.prior.*
    sigma.prior.reference <- sigma.prior.reference[params.id[(1:n)]]
    sigma.prior.probability <- sigma.prior.probability[params.id[(1:n)]]
    sigma.fixed <- sigma.fixed[params.id[(1:n)]]
    nUnkSigmas <- length(sigma.prior.reference)

    if(missing(low.params.fixed)) {
      low.params.fixed <- rep(NA, nEdges)
    } else {
      stopifnot(length(low.params.fixed)==nEdges)
    }
    low.params.fixed[params.id[n+1:nEdges]-n]
    if(any(!is.na(low.params.fixed)))  stop("WORK IN PROGRESS!")

    ii <- c(1:n, qij$ii)
    jj <- c(1:n, qij$jj)
    ii <- ii[order(jj)]
    jj <- jj[order(jj)]

    iuq <- qij$ilq  ## mem order
    iuqpac <- qij$ilqpac

    ifi <- row(Q0)[qij$ifil]
    jfi <- col(Q0)[qij$ifil]

    if(nEdges==0) {
      stop("This graph is trivial, please consider 'iid' model!")
    }
    if(missing(base)){
      warning("Missing base model! Using 'iid'.")
      base <- rep(0, nEdges)
    }

    Ibase <- hessian(model, base)
    if(debug) {
      cat("I(base model) elements\n")
      print(str(Ibase))
    }
    stopifnot(all(dim(Ibase) == c(nEdges, nEdges)))
    ## this is I(\theta_0)^{-0.5} * \theta_0
    thetabasescaled <- drop(attr(Ibase, "hneg.5") %*%
                              attr(Ibase, "base"))

    ## constant
    lc <- log(lambda) + lgamma(1.0+nEdges/2)
    lc <- lc -log(nEdges)-(0.5*nEdges)*log(pi)
    lc <- lc #+ sum(log(attr(Ibase, "decomposition")$values))

    if(debug) {
      cat('log C', lc, '\n')
    }

    cgmodel <- "inla_cgeneric_graphpcor"
    the_model <- list(
      f = list(
        model = "cgeneric",
        n = as.integer(n),
        cgeneric = list(
          model = cgmodel,
          shlib = libpath,
          n = as.integer(n),
          debug = as.integer(debug)),
        data = list(
          ints = list(
            n = as.integer(n),
            debug = as.integer(debug),
            ne = as.integer(nEdges),
            nfi = as.integer(nfi),
            ii = as.integer(jj-1),
            jj = as.integer(ii-1),
            iuq = as.integer(iuq-1),
            iuqpac = as.integer(iuqpac-1),
            ifi = as.integer(ifi-1),
            jfi = as.integer(jfi-1),
            itheta = as.integer(params.id -1),
            sfixed = as.integer(sigma.fixed)),
          doubles = list(
            lambda = as.numeric(lambda),
            sigmaref = as.numeric(sigma.prior.reference),
            sigmaprob = as.numeric(sigma.prior.probability),
            lconst = as.numeric(lc),
            thetabasescaled = as.numeric(thetabasescaled),
            thetab = as.numeric(base)),
          characters = list(
            model = cgmodel,
            shlib = libpath
          ),
          matrices = list(
            hHneg = attr(Ibase, "hneg.5"),
            H = matrix(new("numeric", Ibase), nrow(Ibase))),
          smatrices = list()
        )
      )
    )

    class(the_model) <- c("cgeneric", ## INLAtools
                          "inla.cgeneric") ## this is needed in INLA::f()
    class(the_model$f$cgeneric) <- class(the_model)
    return(the_model)

  }
