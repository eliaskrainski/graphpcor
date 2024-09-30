#' Define cgeneric method for generic0 model
#' @param model an string equal 'generic0'
#' Implement a generic model where the precision matrix is
#'  \deqn{Q = \tau R}
#' where the structure matrix R is supplied by the user
#' and \eqn{\tau} is the (local, see detais) precision parameter.
#' This uses the cgeneric interface that can be used as a
#' model in a `INLA` `f()` model component.
#' @param R the structure matrix for the model definition.
#' @param param length two vector with the parameters
#' \eqn{a,b} for the PC-prior distribution defined from
#'   \deqn{P(\sigma > a) = u}
#' where \eqn{\sigma} can be interpreted as marginal standard
#' deviation of the process if scale = TRUE. See details.
#' @param scale logical indicating if it is to scale
#' the R structure matrix so that the geometric mean for the
#' marginal variances is equal to one when the precision is 1.
#' @param debug logical indicating if it is to debug.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. It is not considered if
#' libpath is provided.
#' @return a [inla.cgeneric] object to be used in the f() formula term in INLA.
#' @details
#' If scale = TRUE the matrix \eqn{R} is scaled so that
#'  \deqn{Q = \tau s R}
#'  where \eqn{s} is the geometric mean of the diagonal
#'  elements of the generalized inverse of \eqn{R}.
#' \deqn{s = \exp{\sum_i \log((R^{-})_{ii})/n}}
#'
cgeneric_generic0 <-
  function(R,
           param,
           constr = TRUE,
           scale = TRUE,
           debug = FALSE,
           useINLAprecomp = !TRUE) {

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

    stopifnot(param[1]>0)
    if(is.na(param[2])) {
      param[2] = 0.0
    }
    stopifnot(param[2]>=0)
    stopifnot(param[2]<=1)

    R <- INLA::inla.as.sparse(R)

    n <- as.integer(nrow(R))
    stopifnot(n>0)

    idx <- which(R@i <= R@j)

    if(debug) {
      print(str(list(
        ii = R@i,
        jj = R@j,
        idx = idx
      )))
    }

    if(scale) {
      R <- INLA::inla.as.sparse(
        inla.scale.model(
          Q = R,
          constr = list(A = matrix(1, 1, n), e = 0)
        )
      )
    }

    ord <- order(R@i[idx])
    nnz <- length(idx)
    cmodel = "inla_cgeneric_generic0"

    the_model <- list(
      f = list(
        model = "cgeneric",
        n = n,
        cgeneric = list(
          model = cmodel,
          shlib = libpath,
          n = n,
          debug = as.logical(debug),
          data = list(
            ints = list(
              n = n,
              debug = debug
            ),
            doubles = list(
              param = param
            ),
            characters = list(
              model = cmodel,
              shlib = libpath
            ),
            matrices = list(
            ),
            smatrices = list(
              Rgraph = c(
                n, n, nnz,
                R@i[idx][ord],
                R@j[idx][ord],
                R@x[idx][ord]
              )
            )
          )
        )
      )
    )

    class(the_model) <- "inla.cgeneric"
    class(the_model$f$cgeneric) <- "inla.cgeneric"

    if(constr) {
      the_model$f$extraconstr <- list(
        A = matrix(1, 1, n),
        e = 0
      )
    }

    return(the_model)

}

#' The 'iid' model uses the 'generic0' with
#' structure matrix as the identity.
#' @param n size of the model
cgeneric_iid <-
  function(n,
           param,
           constr = TRUE,
           scale = TRUE,
           debug = FALSE,
           useINLAprecomp = !TRUE) {
    cgeneric_generic0(
      R = Diagonal(n = n,
                   x = rep(1, n)),
      param = param,
      constr = constr,
      scale = FALSE,
      debug = debug,
      useINLAprecomp = useINLAprecomp
    )
}

