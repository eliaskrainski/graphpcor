#' Build an `cgeneric` to implement the Wishart
#' prior for a precision matrix.
#' @param n integer to define the size of the precision matrix
#' @param dof degrees of freedom model parameter
#' @param R lower triangle of the scale matrix parameter
#' @param ... additional arguments passed on to
#' [INLAtools::cgeneric()], such as `debug`,
#' `shlib` and `useINLAprecomp`.
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
           ...) {

    dotArgs <- list(...)
    if(is.null(dotArgs$debug)) {
      dotArgs$debug <- FALSE
    }

    stopifnot(n>=1)
    stopifnot(dof>(n+1))
    stopifnot(length(R) == (n*(n+1)/2))
    M <- as.integer(length(R))

    if(dotArgs$debug) {
      cat("N = ", n, ", M = ", M, "\n")
    }

    if(is.null(dotArgs$shlib)) {
      if(dotArgs$debug){
        cat("searching shlib...\n")
      }
      dotArgs$shlib <-
        INLAtools::cgeneric_shlib_path(
          package = "graphpcor",
          useINLAprecomp = dotArgs$useINLAprecomp,
          debug = dotArgs$debug
        )
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

    if(dotArgs$debug) {
      cat('hldet = ', hldetr, ', log const = ', lcprior, '\n')
    }

    cmodel = "inla_cgeneric_Wishart"

    the_model <- do.call(
      what = INLAtools::cgenericBuilder,
      args = list(
          model = cmodel,
          n = as.integer(n),
          debug = as.integer(dotArgs$debug), ## 1
          shlib = dotArgs$shlib,
          dof = as.numeric(dof),
          lcprior = as.double(lcprior),
          R = as.numeric(rr)
      ))

    return(the_model)

}
