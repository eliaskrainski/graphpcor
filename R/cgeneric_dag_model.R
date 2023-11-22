#' Build the objects to implement the model using the
#' cgeneric method in `INLA` from a list of expressions
#' defining a Direct Acyclic Graph - DAG correlation model
#' to be used as a model in a `INLA` `f()` model component.
#'
#' @param dag the DAG model.
#' @param lambda the lambda for the graph correlation prior.
#' @param sigma.prior.reference a vector with the reference values
#' to define the prior for the standard deviation parameters.
#' @param sigma.prior.probability a vector with the probability values
#' to define the prior for the standard deviation parameters.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. It is not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @details
#'  The correlation prior as in the paper depends on the lambda value.
#'  The prior for each \eqn{sigma_i} is the Penalized-complexity prior
#' which can be defined from the following probability statement
#'  P(sigma > U) = a.
#' where "U" is a reference value and "a" is a probability.
#' The values "U" and probabilities "a" for each \eqn{sigma_i}
#' are passed in the `sigma.prior.reference` and `sigma.prior.probability`
#' arguments.
#' If a=0 then U is taken to be the fixed value of the corresponding sigma.
#' E.g. if there are three sigmas in the model and one supply
#'  sigma.prior.reference = c(1, 2, 3) and
#'  sigma.prior.probability = c(0.05, 0.0, 0.01)
#' then sigma[2] is fixed to 2 and not estimated.
#' @return objects to be used in the f() formula term in INLA.
#' @export
cgeneric0_dag_model <-
  function(dag,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           debug = FALSE,
           useINLAprecomp = !TRUE,
           libpath = NULL) {

    if (is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("corGraphs")
      } else {
        libpath <- system.file("libs", package = "corGraphs")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "corGraphs.dll")
        } else {
          libpath <- file.path(libpath, "corGraphs.so")
        }
      }
    }

    d.el <- dag_elements(dag)
    q.el <- dag_e2precision(dag)
    NC <- q.el$nc

    stopifnot(length(sigma.prior.reference) == NC)
    stopifnot(length(sigma.prior.probability) == NC)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    qc <- q.el$q[1:NC, 1:NC]
    diag(qc) <- diag(qc) + 1e-9
    ii <- col(qc)[!upper.tri(qc)]
    jj <- row(qc)[!upper.tri(qc)]

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs0",
        shlib = libpath,
        n = as.integer(q.el$nc),
        debug = as.integer(debug),
        np = as.integer(q.el$np),
        i2th = as.integer(q.el$i2th-1L),
        iq2th = as.integer(q.el$iq2th-1L),
        i1th = as.integer(q.el$i1th-1L),
        iq1th = as.integer(q.el$iq1th-1L),
        iq1ch = as.integer(q.el$iq1ch-1L),
        ii = as.integer(ii - 1L),
        jj = as.integer(jj - 1L),
        sch = as.double(q.el$sch),
        sth = as.double(q.el$sth),
        q = as.double(q.el$q),
        lambda = as.double(lambda),
        slambdas = as.double(slambdas)
      )
    )

    return(the_model)
  }
#' Implement the model from covariance
#' @return objects to be used in the f() formula term in INLA.
#' @export
cgeneric_dag_model <-
  function(dag,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           debug = FALSE,
           useINLAprecomp = !TRUE,
           libpath = NULL) {

    if (is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("corGraphs")
      } else {
        libpath <- system.file("libs", package = "corGraphs")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "corGraphs.dll")
        } else {
          libpath <- file.path(libpath, "corGraphs.so")
        }
      }
    }

    d.el <- dag_elements(dag)
    ich <- unlist(lapply(d.el, function(x)
      x$id[!x$parent]))
    sch <- unlist(lapply(d.el, function(x)
      x$signal[!x$parent]))
    sch <- sch[ich]
    d.elc <- dag_e2covariance(d.el)
    np <- length(dag)
    nv <- sapply(d.elc$iv, length)
    iiv <- rep(1:np, nv)
    jjv <- unlist(d.elc$iv)
    itop <- d.elc$itop
    nc <- nrow(itop)
    ii <- col(itop)[!upper.tri(itop)]
    jj <- row(itop)[!upper.tri(itop)]

    stopifnot(length(sigma.prior.reference) == nc)
    stopifnot(length(sigma.prior.probability) == nc)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs_sfixed",
        shlib = libpath,
        n = as.integer(nc),
        debug = as.integer(debug),
        np = as.integer(np),
        nv = as.integer(nv),
        ipar = as.integer(d.elc$iparent-1L),
        iiv = as.integer(iiv-1L),
        jjv = as.integer(jjv-1L),
        itop = as.integer(itop-1L),
        ii = as.integer(ii-1L),
        jj = as.integer(jj-1L),
        lambda = as.double(lambda),
        slambdas = as.double(slambdas),
        schildren = as.double(sch)
      )
    )

    return(the_model)
  }
