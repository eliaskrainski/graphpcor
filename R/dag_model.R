#' Build the objects to implement a
#' Direct Acyclic Graph - DCG correlation model
#' to be used as a model in a `INLA` `f()` model component.
#' @return objects to be used in the f() formula term in INLA.
#' @export
dag_model <-
  function(graph,
           sigma.prior.reference,
           sigma.prior.probability,
           lambda,
           debug = FALSE,
           useINLAprecomp = !TRUE,
           libpath = NULL) {

    return(graph_qchol_index(graph))

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

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs_cholQ",
        shlib = libpath,
        n = as.integer(nc),
        debug = as.integer(debug),
        np = as.integer(np),
        i
      )
    )

    return(the_model)

}
