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

    if(is.matrix(graph))
      graph <- matrix2graph(graph)
    qij <- graph_qchol_index(graph)
    n <- as.integer(qij$n)
    stopifnot(n>0)

    stopifnot(all(lambda>0))
    stopifnot(length(sigma.prior.reference) == n)
    stopifnot(length(sigma.prior.probability) == n)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    ne <- as.integer(length(qij$ii))
    nnz <- n + ne
    nfi <- as.integer(length(qij$ifil))

    ii <- as.integer(c(1:n, qij$ii)-1L)
    jj <- as.integer(c(1:n, qij$jj)-1L)
    ii <- ii[order(jj)]
    jj <- jj[order(jj)]
    ilq <- as.integer(qij$ilq-1L)
    iuq <- as.integer(qij$iuq-1L)

    ifi <- as.integer(row(qij$Q)[qij$ifil]-1L)
    jfi <- as.integer(col(qij$Q)[qij$ifil]-1L)

    the_model <- do.call(
      "inla.cgeneric.define",
      list(
        model = "inla_cgeneric_corgraphs_cholQ",
        shlib = libpath,
        n = n,
        debug = as.integer(debug),
        ne = ne,
        nnz = nnz,
        nfi = nfi,
        ii = jj,
        jj = ii,
        ilq = ilq,
        iuq = iuq,
        ifi = ifi,
        jfi = jfi,
        lambda = lambda,
        slambdas = slambdas
      )
    )

    return(the_model)

}
