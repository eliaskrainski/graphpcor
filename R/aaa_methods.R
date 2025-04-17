#' @rdname graphpcor
#' @title The `graphpcor` generic method for [graphpcor-class]
#' @param ... either a list of formulae or a matrix
#' @return a `graphpcor` object
#' @export
graphpcor <- function(...) {
  UseMethod("graphpcor")
}
#' The `prec` method
#' @rdname prec-methods
#' @param model a model object
#' @param ... additional arguments
#' @return a precision matrix
#' @export
prec <- function(model, ...) {
  UseMethod("prec")
}
#' @describeIn prec-methods
#' The default precision method
#' computes the inverse of the variance
#' @export
prec.default <- function(model, ...) {
  v <- vcov(model, ...)
  return(
    forwardsolve(
      backsolve(
        chol(v)
      )
    )
  )
}
#' @describeIn prec-methods
#' Define the prec method for an inla output object
#' @export
prec.inla <- function(model, ...) {
  if(is.null(model$misc$config$config)) {
    warning("inla.rerun() with config = TRUE in control.compute.")
    model$.args$control.compute$config <- TRUE
    model <- do.call("inla", args = model$.args)
  }
  Qu <- INLA::inla.as.sparse(
    model$misc$config$config[[1]]$Qprior
  )
  #  ii <- which(Qu@i < Qu@j)
  # if(length(ii)>0) {
  Q <- #inla.as.sparse(
    Matrix::sparseMatrix(
      #        i = c(Qu@i, Qu@j[ii]) + 1L,
      #       j = c(Qu@j, Qu@i[ii]) + 1L,
      #      x = c(Qu@x, Qu@x[ii])
      i = Qu@i + 1L,
      j = Qu@j + 1L,
      x = Qu@x,
      symmetric = TRUE,
      repr = "T"
    )
  #    )
  # } else {
  #  Q <- Qu
  #}
  return(Q)
}
#' Define the is.zero method
#' @param x an R object
#' @param ... additional arguments
#' @return logical
#' @export
is.zero <- function(x, ...) {
  UseMethod("is.zero")
}
#' @describeIn is.zero
#' The is.zero.default definition
#' @export
is.zero.default <- function(x, ...) {
  a <- abs(as.numeric(c(x)))
  if(diff(range(a))<(.Machine$double.eps^0.9)) {
    tol <- (.Machine$double.eps^0.9)
  } else {
    tol <- .Machine$double.eps *
      max(sqrt(length(a))) * max(a)
  }
  return(a < tol)
}
#' @describeIn is.zero
#' The is.zero.matrix definition
#' @export
is.zero.matrix <- function(x, ...) {
  stopifnot(inherits(x, "matrix"))
  a <- abs(x)
  if(diff(range(a))<(.Machine$double.eps^0.9)) {
    tol <- (.Machine$double.eps^0.9)
  } else {
    tol <- .Machine$double.eps *
      max(sqrt(length(a))) * max(a)
  }
  return(a < tol)
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
