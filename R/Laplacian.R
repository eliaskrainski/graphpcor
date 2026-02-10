  #' The Laplacian of a graph
#' @rdname Laplacian
#' @param x object defining a graph
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
Laplacian <- function(x) {
  UseMethod("Laplacian")
}
#' @describeIn Laplacian
#' The Laplacian default method (none)
#' @export
Laplacian.default <- function(x) {
  stop("No Laplacian method for", class(x), "!")
}
#' @describeIn Laplacian
#' The Laplacian of a matrix
#' @export
Laplacian.matrix <- function(x) {
  A <- 1 - is.zero(x)
  if(any(A!=t(A)))
    warning("Not symmetric!")
  L <- diag(colSums(A)) - A
  return(L)
}
#' @describeIn Laplacian
#' The Laplacian of a Matrix
#' @importFrom Matrix colSums Diagonal
#' @export
Laplacian.Matrix <- function(x) {
  o <- Sparse(x)
  o@x <- rep(1, length(o@x))
  iid <- which(o@i==o@j)
  if(length(iid)>0) {
    o@x[iid] <- 0
  }
  d <- colSums(o)
  return(Sparse(
    Diagonal(n = nrow(o), x = d)-o
  ))
}
#' @describeIn graphpcor
#' The Laplacian method for a `graphpcor`
#' @export
Laplacian.graphpcor <- function(x) {
  L <- -attr(x, "graph")
  diag(L) <- -colSums(L)
  return(L)
}
