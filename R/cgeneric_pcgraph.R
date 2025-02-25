#' Define the cgeneric model with a PC prior for a
#' Direct Acyclic Graph - DAG correlation model
#' to be used as a model in a `INLA` `f()` model component.
#' @param graph the model graph to define the model structure.
#' @param lambda the parameter for the exponential prior on
#' the radius of the sphere, see details.
#' @param base numeric vector or matrix with the reference
#' correlation model against what the KLD will be evaluated.
#' If it is a matrix, is used as is.
#' If it is a vector, a correlation matrix is defined
#' considering the graph model and this vector as
#' the parameters in the lower triangle matrix L.
#' @param sigma.prior.reference numeric vector with the reference
#' standard deviation to define the PC priors for each one of
#' the marginal variance parameters.
#' @param sigma.prior.probability numeric vector with a probability p
#' to define the PC priors for each one of
#' the marginal variance parameters.
#' #' to make the prior for each marginal standard deviation using
#' P(sigma < `sigma.prior.reference') = p
#' @param debug logical indicating if it is to debug.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. It is not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @return objects to be used in the f() formula term in INLA.
#' @export
cgeneric_pcgraph <-
  function(graph,
           lambda,
           base,
           sigma.prior.reference,
           sigma.prior.probability,
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

    Q0 <- Laplacian(graph)
    n <- nrow(Q0)
    if(debug>99) {
      cat("the 'n ='", n, "graph Laplacian is\n")
      print(Q0)
    }
    l1 <- t(chol(Q0 + diag(1.0, n, n)))
    qnz <- !is.zero(Q0)
    qij <- list(
      ii = row(Q0)[qnz & lower.tri(Q0, diag = FALSE)],
      jj = col(Q0)[qnz & lower.tri(Q0, diag = FALSE)],
      iq = which(Q0!=0))
    qij$ilq <- which(qnz & lower.tri(Q0, diag = TRUE))
    qij$iuq <- which(qnz & upper.tri(Q0, diag = TRUE))
    ll <- t(chol(Q0 + diag(n)))
    qij$ifil <- setdiff(which(ll!=0), qij$ilq)
    if(debug>99) {
      print(qij)
    }
    n <- nrow(Q0)
    stopifnot(n>0)

    stopifnot(all(lambda>0))
    stopifnot(length(sigma.prior.reference) == n)
    stopifnot(length(sigma.prior.probability) == n)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    nEdges <- length(qij$ii)
    nnz <- n + nEdges
    nfi <- length(qij$ifil)

    ii <- c(1:n, qij$ii)
    jj <- c(1:n, qij$jj)
    ii <- ii[order(jj)]
    jj <- jj[order(jj)]

    iuq <- qij$ilq  ## mem order
    ilq <- matrix(1:(n*n), n, n, byrow = TRUE)[qij$ilq]

    ifi <- row(Q0)[qij$ifil]
    jfi <- col(Q0)[qij$ifil]

    if(nEdges==0) {
      stop("This graph is trivial, please consider 'iid' model!")
    }
    if(missing(base)){
      warning("Missing base model! Using 'iid'.")
      base <- rep(0, nEdges)
    }

    H.el <- graph2H(graph, base, method = "eigen")
    if(debug) {
      cat("H elements\n")
      print(str(H.el))
    }

    ## constant
    lc <- log(lambda) -(nEdges-1)*log(pi) - log(2)
    lc <- lc - sum(log(H.el$Hd$values))

    if(debug) {
      cat('log C', lc, '\n')
    }

    m_args <- list(
      model = "inla_cgeneric_pcgraph",
      shlib = libpath,
      n = as.integer(n),
      debug = as.integer(debug),
      ne = as.integer(nEdges),
      nnz = as.integer(nnz),
      nfi = as.integer(nfi),
      ii = as.integer(jj-1),
      jj = as.integer(ii-1),
      ilq = as.integer(ilq-1),
      iuq = as.integer(iuq-1),
      ifi = as.integer(ifi-1),
      jfi = as.integer(jfi-1),
      lambda = as.numeric(lambda),
      slambdas = as.numeric(slambdas),
      lconst = as.numeric(lc)
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
