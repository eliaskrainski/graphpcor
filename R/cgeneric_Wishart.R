#' Build an `cgeneric` to implement the Wishart
#' prior for a precision matrix.
#' @param n integer to define the size of the precision matrix
#' @param dof degrees of freedom model parameter
#' @param R lower triangle of the scale matrix parameter
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'shlib' is provided.
#' @param shlib string, default is NULL, with the path to the shared object.
#' @details
#' For a random \eqn{p\times p} precision matrix \eqn{Q},
#' given the parameters \eqn{d} and \eqn{R},
#' respectively scalar degree of freedom and the _inverse_
#' scale \eqn{p\times p} matrix the Wishart density is
#' \deqn{|Q|^{(d-p-1)/2}\textrm{e}^{-tr(RQ)/2}|R|^{p/2}2^{-dp/2}\Gamma_p(n/2)^{-1}}
#'
#' @return a `cgeneric`, [cgeneric()] object.
cgeneric_Wishart <-
  function(n,
           dof,
           R,
           debug = FALSE,
           useINLAprecomp = TRUE,
           shlib = NULL) {

    if(is.null(shlib)) {
      if (useINLAprecomp) {
        shlib <- INLA::inla.external.lib("graphpcor")
      } else {
        shlib <- system.file("libs", package = "graphpcor")
        if (Sys.info()["sysname"] == "Windows") {
          shlib <- file.path(shlib, "graphpcor.dll")
        } else {
          shlib <- file.path(shlib, "graphpcor.so")
        }
      }
    }
    stopifnot(file.exists(shlib))

    stopifnot(n>=1)
    stopifnot(dof>(n+1))
    stopifnot(length(R) == (n*(n+1)/2))
    M <- as.integer(length(R))

    if(debug) {
      cat("N = ", n, ", M = ", M, "\n")
    }

    rr <- diag(R[1:n], nrow = n, ncol = n)
    if(n>1) {
      il <- (n+1):M
      rr[lower.tri(rr, diag = FALSE)] <- R[il]
      rr <- t(rr)
      rr[lower.tri(rr, diag = FALSE)] <- R[il]
      ur <- chol(rr)
      hldetr <- sum(log(diag(ur)))
    } else {
      hldetr <- log(rr[1])
    }

    lcprior <- dof*hldetr -
      0.5*dof*n*log(2.0) -
      0.25*n*(n-1)*log(pi) -
      sum(lgamma((dof+1-(1:n))/2))

    if(debug) {
      cat('hldet = ', hldetr, ', log const = ', lcprior, '\n')
    }

    cmodel = "inla_cgeneric_Wishart"

    the_model <- do.call(
      what = INLAtools::cgenericBuilder,
      args = list(
        model = cmodel,
        n = as.integer(n),
        debug = as.logical(debug),
        shlib = shlib,
        dof = as.numeric(dof),
        lcprior = as.double(lcprior),
        R = as.numeric(rr)
      ))

    return(the_model)

}
