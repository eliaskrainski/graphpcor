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
           useINLAprecomp = TRUE,
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

    nS <- length(S)
    stilde <- sapply(S, function(x)
      strsplit(as.character(x), split = "~"))
    if(debug)
      print(stilde)
    stopifnot(all(substr(stilde[2, ], 1, 1) == "p"))
    iParents1 <- as.integer(substring(
      unlist(stilde[2, ]), 2))
    if(debug)
      print(iParents1)
    NP <- length(iParents1)
    if(debug)
      cat("NP = ", NP, "\n")
    elements2 <- lapply(strsplit(
      gsub(" ", "", unlist(stilde[3, ])),
      "+", fixed = TRUE), unique)
    if(debug)
      print(elements2)
    stopifnot(all(unique(
      unlist(lapply(elements2, substr, 1, 1))
      ) %in% c("p", "c")))
    iielements2 <- lapply(elements2, function(x)
      as.integer(substring(x, 2)))
    if(debug)
      print(iielements2)
    iparent2 <- lapply(elements2, function(x)
      substr(x, 1, 1) == "p")
    if(debug)
      print(iparent2)
    stopifnot(all(unlist(
      iielements2)[unlist(iparent2)]) %in% iParents1)
    ichildren2 <- lapply(elements2, function(x)
      substr(x, 1, 1) == "c")
    if(debug)
      print(ichildren2)
    NC <- sum(unlist(ichildren2))
    if(debug)
      cat("NC =", NC, "\n")
    stopifnot(max(unlist(
      iielements2)[unlist(ichildren2)]) == NC)
    iiminus <- gregexpr("-", stilde[3,], fixed = TRUE)
    if(debug)
      print(iiminus)
    iiplus <- gregexpr("+", stilde[3,], fixed = TRUE)
    if(debug)
      print(iiplus)
    iisignal <- vector('list', nS)
    for(k in 1:nS) {
      ss <- c(sum(iiminus[[k]]>0),
              sum(iiplus[[k]]>0))
      a <- integer(max(iiminus[[k]], iiplus[[k]]))
      s <- integer()

    }

    if(FALSE) {

      q2i <- c(.1, 100, 10000)
      q2i

      Q <- cbind(rbind(diag(NC), matrix(0, NP, NC)),
               rbind(matrix(0, NC, NP), diag(q2i)))
      Q
      jj <- NC + iParents1
      for(k in 1:nS) {
        if (any(ichildren2[[k]])) {
          ii <- iielements2[[k]][ichildren2[[k]]]
          cat("k1 =", k, "ii =", ii, "\n")
          Q[jj[k], jj[k]] <- Q[jj[k], jj[k]] + length(ii)
          for(i in ii) {
            if(iiminus[[k]])
            Q[i, jj[k]] <- -1
          }
        }
        if (any(iparent2[[k]])) {
          i1 <- NC + iParents1[k]
          ii <- NC + iielements2[[k]][iparent2[[k]]]
          cat("k2 =", k, "ii =", ii, "\n")
          for(i in ii) {
            Q[i1, i1] <- Q[i1, i1] + q2i[i-NC]
            Q[i, i1] <- Q[i, i1] - q2i[i-NC]
          }
        }
      }
      Q

    }

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs",
        shlib = libpath,
        n = n,
        debug = as.integer(debug),
        verbose = as.integer(verbose),
        Rmanifold = as.integer(Rmanifold),
        dimension = as.integer(dimension),
        aaa = as.integer(alphas),
        nm = as.integer(nm),
        cc = as.double(cc),
        bb = mm$bb,
        prs = control.priors$prs,
        prt = control.priors$prt,
        psigma = control.priors$psigma,
        ii = lmats$graph@i,
        jj = lmats$graph@j,
        tt = t(mm$TT),
        xx = t(lmats$xx)
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
