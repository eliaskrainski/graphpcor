#' Build an `cgeneric` for a graph, see [graphpcor()]
#' @description
#' From either a `graphpcor` (see [graphpcor()]) or
#' a square matrix (used as a graph),
#' creates an `cgeneric` (see [INLAtools::cgeneric()])
#' to implement the Penalized Complexity prior using the
#' Kullback-Leibler divergence - KLD from a base graphpcor.
#' @param model a `graphpcor` (see [graphpcor()]) or
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
#' @param sigma.prior.reference numeric vector to set the reference
#' for each standard deviation parameter for its PC-prior.
#' If missing, the model will be assumed as for a correlation.
#' NOTE: `iparams` will be applied here as
#' `sigma.prior.reference[iparams[1:n]]`.
#' @param sigma.prior.probability numeric vector with to
#'set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference`) = p. If missing, all the
#' marginal variances are considered as known.
#' If a vector is given and a probability is NA, 0 or 1, the
#' corresponding `sigma.prior.reference` will be used as fixed.
#' NOTE: `iparams` will be applied here as
#' `sigma.prior.probability[iparams[1:n]]`.
#' @param iparams integer vector with length equal `n+m`,
#' where `m` is the number of correlation parameters,
#' to identify (possible) common parameters in the model
#' @param ifixed integer vector to specify which correlation
#' parameters are treated as known when computing the Hessian.
#' By default all correlation parameters are treated as unknown.
#' Note: it will be applied to `iparams` as `iparams[n+1:m]`.
#' Example: if `ifixed = c(1,3)`, the first and third
#' correlation parameters will be treated as fixed and the
#' Hessian will be computed for the second correlation parameter.
#' Please see the examples in [basepcor()].
#' @param ... additional arguments passed on to
#' [INLAtools::cgeneric()], such as `debug`,
#' `shlib` and `useINLAprecomp`.
#' @seealso [graphpcor()] and [basepcor()]
#' @returns `cgeneric` object.
cgeneric_graphpcor <-
  function(model,
           lambda,
           base,
           sigma.prior.reference,
           sigma.prior.probability,
           iparams,
           ifixed,
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

    if(inherits(model, "matrix")) {
      if(dotArgs$debug) {
        cat("Building 'graphpcor' from a 'matrix'!")
      }
      model <- graphpcor(model)
    }
    if(inherits(model, "Matrix")) {
      if(dotArgs$debug) {
        cat("Building 'graphpcor' from a 'Matrix'!")
      }
      model <- graphpcor(model)
    }

    n <- dim(model)[1]
    stopifnot(n>0)
    nEdges <- dim(model)[2]
    if(nEdges==0) {
      stop("This graph is trivial, please consider 'iid' model!")
    }

    ## Matrix structure
    Q0 <- Laplacian(model)
    if(dotArgs$debug) {
      print(model)
      cat("Laplacian is\n")
      print(Q0)
    }

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

    nnz <- n + nEdges
    nfi <- length(qij$ifil)

    ii <- c(1:n, qij$ii)
    jj <- c(1:n, qij$jj)
    ii <- ii[order(jj)]
    jj <- jj[order(jj)]

    iuq <- qij$ilq  ## mem order
    iuqpac <- qij$ilqpac

    ifi <- row(Q0)[qij$ifil]
    jfi <- col(Q0)[qij$ifil]

    iparams <- m_iparams_fncheck(nnz, iparams)
    npars1 <- length(unique(iparams[1:n]))
    npars2 <- length(unique(iparams[n+1:nEdges]))

    if(dotArgs$debug) {
      print(list(n = n, nnz = nnz, nfi = nfi,
                 itheta = itheta, iparams = iparams,
                 npars1 = npars1, npars2 = npars2))
    }

    ## sigma prior
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
    if(dotArgs$debug) {
      print(list(sigmaref = sigma.prior.reference,
                 sigmaprob = sigma.prior.probability,
                 sfixed = sigma.fixed,
                 nUnkSigmas = nUnkSigmas))
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

    if(missing(base)){
      warning("Missing base model! Using 'iid'.")
      base <- rep(0, length(itheta))
    }

    basemodel <- basepcor(
      base,
      p = n,
      itheta = itheta,
      d0 = n:1,
      iparams = iparams[n+1:nEdges]-npars1
    )
    stopifnot(length(basemodel$theta) == npars2)
    if(dotArgs$debug) {
      print(basemodel)
    }

    theta0 <- basemodel$theta
    I0 <- hessian(basemodel)
    if(dotArgs$debug) {
      cat("Hessian\n")
      print(I0)
    }

    if(is.null(I0)) {
      I0d <- list(logDeterminant = 0,
                  sqrt = matrix(0, 0))
    } else {
      I0d <- dspd(I0)
    }
    if(dotArgs$debug) {
      cat("Hessian decomposition\n")
      print(I0d)
    }

    if(missing(ifixed)) {
      ifixed <- NULL
    }
    cfixed <- ((1:npars2) %in% ifixed)

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
        cfixed = as.integer(cfixed),
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
