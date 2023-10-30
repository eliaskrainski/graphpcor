#' Define a Graph based correlation graph model object
#'  to define a random model for the `INLA` `f()` call.
#'
#' @param graph a graph model.
#' @param lambda the lambda for the graph prior.
#' @param sigma.prior.reference a vector with the reference values
#' to define the prior for the standard deviation parameters.
#' @param sigma.prior.probability a vector with the probability values
#' to define the prior for the standard deviation parameters.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. Not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @details
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
CorGraphs_cgeneric <-
  function(graph,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    if (is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("INLAcorrel")
      } else {
        libpath <- system.file("libs", package = "INLAcorrel")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "INLAcorrel.dll")
        } else {
          libpath <- file.path(libpath, "INLAcorrel.so")
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
          Q[ii, jj[k]] <- -1
          Q[jj[k], jj[k]] <- Q[jj[k], jj[k]] + length(ii)
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
        model = "inla_cgeneric_sstspde",
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
