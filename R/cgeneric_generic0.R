#' Implement the Besag model using cgeneric interface
#' to be used as a model in a `INLA` `f()` model component.
#' @param graph the graph for the model definition.
#' @param param the parameters for the PC-prior distribution
#' on the precision parameter stated as
#'   P(sigma > param[1]) = param[2].
#' @param scale logical indicating if it is to scale
#' the R structure matrix so that the geometric mean for the
#' marginal variances is equal to one when the precision is 1.
#' @param debug logical indicating if it is to debug.
#' @param useINLAprecomp logical indicating if is to be used
#' shared object pre-compiled by INLA. It is not considered if
#' libpath is provided.
#' @param libpath string to the shared object. Default is NULL.
#' @return objects to be used in the f() formula term in INLA.
#' @export
cgeneric_besag <-
  function(graph,
           param,
           constr = TRUE,
           scale = TRUE,
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

    graph <- INLA::inla.read.graph(graph)
    n <- as.integer(graph$n)
    stopifnot(n>0)

    stopifnot(param[1]>0)
    if(is.na(param[2])) {
      param[2] = 0.0
    }
    stopifnot(param[2]>=0)
    stopifnot(param[2]<=1)

    ii <- rep(1:n, graph$nnbs)
    jj <- unlist(graph$nbs[graph$nnbs>0])
    ijs <- which(ii < jj)

    if(debug) {
      print(str(list(
        ii = ii,
        jj = jj,
        ijs = ijs
      )))
    }

    ii <- c(1:n, ii[ijs])
    jj <- c(1:n, jj[ijs])

    R <- INLA::inla.as.sparse(
      sparseMatrix(
        i = ii,
        j = jj,
        x = c(graph$nnbs,
              rep(-1, length(ii)-n)
              ),
        dims = c(n, n)
        )
      )

    if(scale) {
      R <- INLA::inla.as.sparse(
        inla.scale.model(
          Q = R,
          constr = list(A = matrix(1, 1, n), e = 0)
        )
      )
    }

    idx <- which(R@i <= R@j)
    ord <- order(R@i[idx])
    nnz <- length(idx)
    cmodel = "inla_cgeneric_besag"

    the_model <- list(
      f = list(
        model = "cgeneric",
        n = n,
        cgeneric = list(
          model = cmodel,
          shlib = libpath,
          n = n,
          debug = as.integer(debug),
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

    return(the_model)

}
