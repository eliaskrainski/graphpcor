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
#' g <- list(c1~c2, c2~c3, c3~c4, c4~c5, c5~c6)
#' dag_L(g) ## random
#' dag_L(g) ## random
#' dag_L(g,
#'   theta = c(rep(log(3), 6), ## diagonal
#'             rep(-1, 5)) ## off-diagonal
#' )
dag_L <- function(g, theta = NULL) {
  loadNamespace("Matrix")
  if(any(inherits(g, c("matrix", "Matrix"), which = TRUE))) {
    R <- as.matrix(g)
    n <- nrow(R)
    not0 <- (R != 0) & lower.tri(R, diag = FALSE)
    ii <- row(R)[not0]
    jj <- col(R)[not0]
  } else {
    R <- as.matrix(graph_Laplacian(g))
    n <- nrow(R)
    ii <- attr(R, "ii")
    jj <- attr(R, "jj")
  }
  stopifnot(n == max(ii,jj))
  nnzq <- length(ii)
  if(is.null(theta))
    theta <- rnorm(n + nnzq)
  L <- diag(exp(theta[1:n]))
  if(nnzq==0)
    return(INLA::inla.as.sparse(L))
  sL <- t(chol(R + diag(n)))
  lfi <- ((sL!=0) & (R==0)) * 1L
  retL <- as.matrix(sparseMatrix(
    i = ii,
    j = jj,
    x = theta[n + 1:nnzq],
    dims = c(n, n)
  )) + L
  if(sum(lfi)>0) {
    retL <- fiL(retL, lfi)
  }
  return(INLA::inla.as.sparse(retL))
}
