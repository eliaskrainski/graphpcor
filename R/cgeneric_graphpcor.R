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
#' to identify (possible) common parameters in the model.
#' Default is `1:(n+m)`.
#' Note: `c(1,2,1)` is allowed, but `c(2,1,2)` is not.
#' @param cfixed integer vector to specify which correlation
#' parameters are treated as known and fixed.
#' By default all correlation parameters are treated as unknown.
#' Example: if `cfixed = c(1,3)`, the first and third
#' correlation parameters will be treated as fixed and the
#' Hessian will be computed for the second correlation parameter.
#' Please see the examples in [basepcor()].
#' Note: consider `iparams[n+1:m]-iparams[n]`.
#' @param d0 numeric vector to specify the diagonal of the
#' Cholesky factor for the initial precision matrix `Q0`,
#'  passed on to [basepcor()]. Default is `d0 = n:1`.
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
           cfixed,
           d0,
           ...) {

    ## collect ... args and check debug
    dotArgs <- list(...)
    if(!any(names(dotArgs)=="debug")) {
      dotArgs$debug <- FALSE
    }

    if(is.null(dotArgs$useINLAprecomp) ||
       dotArgs$useINLAprecomp) {
      vs <- "26.02.06" ## after code structure update
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
      dotArgs$shlib <- do.call(
        what = INLAtools::cgeneric_shlib,
        args = list(debug = dotArgs$debug,
                    package = "graphpcor",
                    useINLAprecomp = dotArgs$useINLAprecomp)
      )
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
      stop("This graph has no edges, consider 'iid' model!")
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

    ## note: iparams for sigma and correl params
    iparams <- m_iparams_fncheck(nnz, iparams)
    npars1 <- length(unique(iparams[1:n]))
    stopifnot(npars1>0)
    npars2 <- length(unique(iparams[n+1:nEdges]))
    stopifnot(npars2>0)

    if(dotArgs$debug) {
      print(list(n = n, nnz = nnz, nfi = nfi,
                 itheta = itheta, iparams = iparams,
                 npars1 = npars1, npars2 = npars2))
    }

    if(missing(base)){
      warning("Missing base model! Using 'iid'.")
      base <- rep(0, npars2)
    }

    if(missing(d0)) {
      d0 <- n:1
    } else {
      if(is.null(dotArgs$useINLAprecomp) ||
         dotArgs$useINLAprecomp) {
        vs <- "26.02.05" ## after code structure update
        ivs <- packageCheck(
          name = "INLA",
          minimum_version = vs) >= vs
        if(is.na(ivs)) {
          warning("'d0' is supported after 'INLA 26-02-03'!")
          stop("Update INLA or use 'useINLAprecomp = TRUE'\n")
        }
      }
    }
    basemodel <- basepcor(
      base,
      p = n,
      itheta = itheta,
      d0 = d0,
      iparams = iparams[n+1:nEdges]-npars1
    )
    stopifnot(length(basemodel$theta) == npars2)
    if(dotArgs$debug) {
      print(basemodel)
    }

    theta0 <- basemodel$theta
    if(missing(cfixed)) {
      cfixed <- NULL
    }
    cfixed <- ((1:npars2) %in% cfixed)
    I0 <- hessian(basemodel, ifixed = which(cfixed))
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

    # ## update sigma.* with iparams (left this for C code)
    # sigma.prior.reference <- sigma.prior.reference[iparams[(1:n)]]
    # sigma.prior.probability <- sigma.prior.probability[iparams[(1:n)]]
    # sigma.fixed <- sigma.fixed[iparams[(1:n)]]
    # if(any(is.na(sigma.fixed)))
    #   stop("Some error in `iparams` or  `sigma.fixed`!")
    # if(dotArgs$debug) {
    #   cat("expanded:\n")
    #   print(list(sigmaref = sigma.prior.reference,
    #              sigmaprob = sigma.prior.probability,
    #              sfixed = sigma.fixed))
    # }

    the_model <- do.call(
      what = INLAtools::cgenericBuilder,
      args = list(
        model = "inla_cgeneric_graphpcor",
        n = as.integer(n),                 ## 0
        debug = as.integer(dotArgs$debug), ## 1
        shlib = dotArgs$shlib,
        ne = as.integer(nEdges),  ## 2
        nfi = as.integer(nfi),    ## 3
        ii = as.integer(jj-1),    ## 4
        jj = as.integer(ii-1),    ## 5
        iuq = as.integer(iuq-1),  ## 6
        iuqpac = as.integer(iuqpac-1), ## 7
        ifi = as.integer(ifi-1),       ## 8
        jfi = as.integer(jfi-1),       ## 9
        ifixed = as.integer(c(sigma.fixed, cfixed)), ## 10
        itheta = as.integer(iparams-1),              ## 11
        lambda = as.numeric(lambda),                      ## 0
        sigmaref = as.numeric(sigma.prior.reference),     ## 1
        sigmaprob = as.numeric(sigma.prior.probability),  ## 2
        lconst = as.numeric(I0d$logDeterminant),          ## 3
        thetabase = as.numeric(theta0),                   ## 4
        d0 = as.numeric(basemodel$d0),                    ## 5
        Ihalf = I0d$sqrt
      )
    )
    return(the_model)

  }
