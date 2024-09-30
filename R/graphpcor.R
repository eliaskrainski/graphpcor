#' Build the objects needed to implement a
#' Direct Cyclic Graph - DCG and
#' Direct Acyclic Graph - DAG correlation models
#' to be used as a model in a `INLA` `f()` model component.
#' @useDynLib corGraphs, .registration = TRUE
#' @importFrom Matrix t
#' @export
corGraphs <- function() {
  cat("Welcome!\n")
}
