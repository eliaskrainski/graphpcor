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
cgeneric_dag_model <-
  function(dag,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
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

    q.el <- dag_precision_elements(dag)
    NC <- q.el$n

    stopifnot(length(sigma.prior.reference) == NC)
    stopifnot(length(sigma.prior.probability) == NC)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs",
        shlib = libpath,
        n = as.integer(q.el$n),
        p = as.integer(q.el$p),
        i2th = as.integer(q.el$i2th),
        iq2th = as.integer(q.el$iq2th),
        i1th = as.integer(q.el$i1th),
        iq1th = as.integer(q.el$iq1th),
        iq1ch = as.integer(q.el$iq1ch),
        sch = as.double(q.el$sch),
        sth = as.double(q.el$sth),
        q = as.double(q.el$q0),
        lambda = as.double(lambda),
        slambdas = as.double(slambdas)
      )
    )

    return(the_model)
  }
