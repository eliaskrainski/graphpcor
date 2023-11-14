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

    el <- dag_elements(dag)
    stopifnot(all(substr(names(el), 1, 1) == "p"))
    ip <- as.integer(substring(names(el), 2))
    stopifnot(length(ip) == length(unique(ip)))
    p <- length(ip)
    n <- sum(sapply(el, function(x) sum(!x$parent)))
    p.nc <- sapply(el, function(x) x$n)
    dd <- c(rep(1, n), p.nc)
    stopifnot((n+p) == length(dd))
    q0 <- diag(x = dd, nrow = n + p, ncol = n + p)
    ij <- matrix(1:((n+p)^2), n+p, n+p)
    iq1th <- integer(2 * (p - 1))
    sth <- i1th <- integer(p-1)
    iq2th <- i2th <- integer(p)
    k2 <- k1 <- 0
    for(i in 1:p) {
      i0 <- which(!el[[i]]$parent)
      if(length(i0)>0) {
        j <- el[[i]]$id[i0]
        q0[j, n+i] <- -el[[i]]$signal[i0]
        q0[n+i, j] <- -el[[i]]$signal[i0]
      }
      i2th[k1 + 1] <- i
      iq2th[k1 + 1] <- ij[(col(ij) == (n+i)) & (row(ij) == (n+i))]
      k1 <- k1 + 1
      i0 <- which(el[[i]]$parent)
      nj <- length(i0)
      if(nj>0) {
        j0 <- el[[i]]$id[i0]
        i1th[k2 + 1:nj] <- j0
        sth[k2 + 1:nj] <- el[[i]]$signal[i0] ## carry on the signal
        j <- n + j0
        iq1th[k2 + 1:nj] <- ij[, n+i][j]
        k2 <- k2 + nj
        i1th[k2 + 1:nj] <- j0
        sth[k2 + 1:nj] <- el[[i]]$signal[i0] ## carry on the signal
        iq1th[k2 + 1:nj] <- ij[n+i, ][j]
        k2 <- k2 + nj
      }
    }
    stopifnot(k1 == p)
    stopifnot(k2 == (2*(p-1)))

    stopifnot(length(sigma.prior.reference) == n)
    stopifnot(length(sigma.prior.probability) == n)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs",
        shlib = libpath,
        n = as.integer(n),
        p = as.integer(p),
        i2th = as.integer(i2th),
        iq2th = as.integer(iq2th),
        i1th = as.integer(i1th),
        sth = as.integer(sth),
        iq1th = as.integer(iq1th),
        q = as.double(q0),
        lambda = lambda,
        slambdas = slambdas
      )
    )
    if (constr) {
      the_model$f$extraconstr <- mm$extraconstr
    }
    # Prepend specialised model class identifier, for bru_mapper use:
    class(the_model) <- c("stModel_cgeneric", class(the_model))
    # Add objects needed by bru_get_mapper.stModel_cgeneric:
    # (alternatively, construct the mapper already here, but that would
    # require loading inlabru even when it's not going to be used)
    the_model[["smesh"]] <- smesh
    the_model[["tmesh"]] <- tmesh

    the_model
  }
