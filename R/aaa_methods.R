#' @rdname graphpcor
#' @title The `graphpcor` generic method for [graphpcor-class]
#' @param ... either a list of formulae or a matrix
#' @return a `graphpcor` object
#' @export
graphpcor <- function(...) {
  UseMethod("graphpcor")
}
#' The Laplacian of a graph
#' @rdname Laplacian
#' @param graph object defining a graph
#' @description
#' The (symmetric) Laplacian of a graph is a
#' square matrix with dimention
#' equal the number of nodes.
#' It is defined as
#' \deqn{L_{ij} = n_i \textrm{ if } i=j, -1 \textrm{ if } i\sim j, 0 \textrm{ otherwise}}{%
#'       Lij = ni if i=j, -1 if i~j or 0 otherwise}
#'  where i~j means that there is an edge
#'  between nodes i and j and
#'  n_i is the number of edges including node i.
#' @return matrix as the Laplacian of a graph
#' @export
Laplacian <- function(graph) {
  UseMethod("Laplacian")
}
#' @describeIn Laplacian
#' The Laplacian default method (none)
#' @export
Laplacian.default <- function(graph) {
  stop("No Laplacian for this graph!")
}
#' @rdname pcor
#' @title Information needed for a prior that penalizes
#' the divergence from a base correlation matrix.
#' @returns a pcor object
#' @export
pcor <- function(base, p, parametrization,
                 itheta, d0) {
  UseMethod("pcor")
}
