#' Build an `inla.cgeneric` for a graph, see [corgraph()]
#' @description
#' From either a `graph` (see [graph()]) or
#' a square matrix (used as a graph),
#' creates an `inla.cgeneric` (see [cgeneric()])
#' to implement the Penalized Complexity prior using the
#' Kullback–Leibler divergence - KLD from a base model.
#' @param model  a `graph` (see [graph()]) or
#' a square matrix (used as a graph)
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
#' reference  standard deviation to define the PC prior for each
#' marginal variance parameters.
#' @param sigma.prior.probability numeric vector with length `n`
#' to set the probability statement of the PC prior for each
#' marginal variance parameters. The probability statement is
#' P(sigma < `sigma.prior.reference') = p.
#' If a probability is NA, 0 or 1, the corresponding
#' `sigma.prior.reference` would be taken as fixed.
#' @param params.id integer vector with length `n+m` to specify
#' common parameter values. The first `n` are index to the
#' standard deviations and the remaining `m`
#' are related to partial correlations.
#' Example: By setting params.id = c(1,1,2,3, 4,5,5,6),
#' the first two standard deviations are common and the
#' second and third partial correlations are common as well,
#' giving 6 unknown parameters in the model.
#' @param low.params.fixed logical vector of length `m`
#' to provide the value at which the parameters in the lower
#' of the L matrix are to be fixed. NA indicates not fixed.
#' Example: with low.params.fixed = c(NA, -1, NA, 1) the first
#' and the third of these parameters will be estimated while
#' the second is fixed and equal to -1 and the forth is fixed
#' and equal to 1. NOTE: `params.id` will be applied here as
#' `low.params.fixed[params.id[n+1:m]-n]`
#' @param debug logical indicating if it is to debug.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. It is not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @return objects to be used in the f() formula term in INLA.
#' @export
cgeneric_corgraph <-
  function(model,
           lambda,
           base,
           sigma.prior.reference,
           sigma.prior.probability,
           params.id,
           low.params.fixed,
           debug = FALSE,
           useINLAprecomp = !TRUE,
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
      model <- corgraph(model)
    }
    Q0 <- Laplacian(model)
    n <- nrow(Q0)
    stopifnot(n>0)
    if(debug>99) {
      cat("the 'n ='", n, "graph Laplacian is\n")
      print(Q0)
    }
    stopifnot(all(lambda>0))
    stopifnot(length(sigma.prior.reference) == n)
    stopifnot(length(sigma.prior.probability) == n)
    sigma.prior.probability[is.na(sigma.prior.probability)] <- 0
    stopifnot(all(sigma.prior.probability>=0.0))
    stopifnot(all(sigma.prior.probability<=1.0))
    sigma.fixed <- is.zero(sigma.prior.probability) |
      is.zero(1-sigma.prior.probability)
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

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
    }
    if(params.id[nnz]<nnz) stop("WORK IN PROGRESS!")

    if(missing(low.params.fixed)) {
      low.params.fixed <- rep(NA, nEdges)
    } else {
      stopfinot(length(low.params.fixed)==nEdges)
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

    Ibase <- hessian(model, base, decomposition = "eigen")
    if(debug) {
      cat("I(base model) elements\n")
      print(str(Ibase))
    }
    stopifnot(all(dim(Ibase) == c(nEdges, nEdges)))
    ## this is I(\theta_0)^{-0.5} * \theta_0
    thetabasescaled <- drop(attr(Ibase, "hneg.5") %*%
                              attr(Ibase, "base"))

    ## constant
    lc <- log(lambda) -(nEdges-1)*log(pi) - log(2)
    lc <- lc - sum(log(attr(Ibase, "decomposition")$values))

    if(debug) {
      cat('log C', lc, '\n')
    }

    m_args <- list(
      model = "inla_cgeneric_corgraph",
      shlib = libpath,
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
      lambda = as.numeric(lambda),
      slambdas = as.numeric(slambdas),
      lconst = as.numeric(lc),
      thetabasescaled = as.numeric(thetabasescaled),
      hHneg = attr(Ibase, "hneg.5")
    )

    if(debug>9) {
      print(str(m_args))
    }

    the_model <- do.call(
      "inla.cgeneric.define",
      m_args
    )

    return(the_model)

  }
