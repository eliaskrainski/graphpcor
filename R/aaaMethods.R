#' Compute the precision matrix
#' @param x object or model
#' @param ... additional arguments passed on
#' for some methods.
#' @export
#' @rdname precision
#' @export
precision <- function(x, ...) {
  UseMethod("precision")
}
#' The default method for precision
#' compute the inverse of the [variance()]
#' @rdname precision
#' @seealso [variance()]
#' @export
precision.default <- function(x, ...) {
  v <- variance(x, ...)
  return(
    forwardsolve(
      backsolve(
        chol(v)
      )
    )
  )
}
#' @describeIn precision
#' Define the precision method for an inla output object
#' @export
precision.inla <- function(x, ...) {
  if(is.null(x$misc$config$config)) {
    warning("inla.rerun() with config = TRUE in control.compute.")
    x$.args$control.compute$config <- TRUE
    x <- do.call("inla", args = x$.args)
  }
  Qu <- INLA::inla.as.sparse(
    x$misc$config$config[[1]]$Qprior
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
#' The variance method
#' @rdname variance
variance <- function(x, ...) {
  UseMethod("variance")
}
#' The variance default method
#' @rdname variance
#' @export
variance.default <- function(x, ...) {
  return(var(x))
}
#' Define the is.zero method
#' @export
is.zero <- function(x, ...) {
  UseMethod("is.zero")
}
#' The is.zero.default definition
#' @export
is.zero.default <- function(x, ...) {
  a <- abs(as.numeric(c(x)))
  tol <- .Machine$double.eps *
    max(sqrt(length(a))) * max(a)
  return(a < tol)
}
#' The is.zero.matrix definition
#' @export
is.zero.matrix <- function(x, ...) {
  stopifnot(inherits(x, "matrix"))
  a <- abs(x)
  tol <- .Machine$double.eps *
    max(dim(a)) * max(a)
  return(a < tol)
}
#' The Laplacian of a graph
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
#' @export
Laplacian <- function(graph) {
  UseMethod("Laplacian")
}
#' The Laplacian.default method
#' @export
Laplacian.default <- function(graph) {
  stop("No Laplacian for this object!")
}
