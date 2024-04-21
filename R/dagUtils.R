#' Function to fill-in Cholesky a matrix
#' @param L the lower triangle of the Cholesky decomposition
#' @param lfi indicator of fill-in elements
#' @return a filled L matrix.
fiL <- function(L, lfi) {
  L <- as.matrix(L)
  ii <- row(lfi)[as.matrix(lfi)!=0]
  jj <- col(lfi)[as.matrix(lfi)!=0]
  for(v in 1:length(ii)) {
    i <- ii[v]
    j <- jj[v]
    if(j==1) {
      warning("j = 1!\n")
      L[i,1] <- 0.0
    } else {
      stopifnot((i>1) & (j>1)) ## L_{11} not allowed
      stopifnot(j>1) ## j=1 is not allowed
      k <- 1:(j-1)
      L[i, j] <- -sum(L[i, k] * L[j, k]) / L[j, j]
    }
  }
  return(L)
}
#' Function to build a precision matrix from a graph model
#' @param g the graph model (as a square matrix)
#' @param theta the vector of parameter in internal scale.
#' Default is NULL and it will be random in this case.
#' @export
#' @examples
#' g <- sparseMatrix(
#'    i = c(1, 2, 3, 4, 5),
#'    j = c(3, 3, 5, 5, 6),
#'    x = 1L,
#'    dims = c(6, 6)
#' ); g <- g + t(g)
#' dag_L(g) ## random
#' dag_L(g) ## random
#' dag_L(g,
#'   theta = c(rep(log(3), 6), ## diagonal
#'             rep(-1, 5)) ## off-diagonal
#' )
dag_L <- function(g, theta = NULL) {
  loadNamespace("Matrix")
  g <- INLA::inla.as.sparse(g)
  n <- nrow(g)
  patt <- as.matrix((g !=0 ) + (Matrix::t(g) != 0))
  patt[patt>1] <- 1L
  R <- INLA::inla.as.sparse(
    Diagonal(n, 1L + rowSums(patt)) -(patt))
  ij <- which(R@i > R@j)
  nnzq <- length(ij)
  if(is.null(theta))
    theta <- rnorm(n + nnzq)
  L <- Diagonal(n, exp(theta[1:n]))
  if(nnzq==0)
    return(crossprod(L))
  sL <- t(chol(R))
  retL <- as.matrix(sparseMatrix(
    i = R@i[ij] + 1L,
    j = R@j[ij] + 1L,
    x = theta[n + 1:nnzq],
    dims = c(n, n)
  ) + L)
  lfi <- ((sL!=0) & (R==0)) * 1L
  if(sum(lfi)>0) {
    retL <- fiL(retL, lfi)
  }
  return(INLA::inla.as.sparse(retL))
}
