#' Function to fill-in Cholesky a matrix
#' @param L the lower triangle of the Cholesky decomposition
#' @param lfi indicator of fill-in elements
#' @return a filled L matrix.
fiL <- function(L, lfi) {
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
#' g <- inla.as.sparse(
#'   sparseMatrix(
#'    i = c(1, 2, 3, 4, 5),
#'    j = c(3, 3, 5, 5, 6),
#'    x = 1L,
#'    dims = c(6, 6)
#'   )
#' )
#' dag_L(g)
dag_L <- function(g, theta = NULL) {
  g <- inla.as.sparse(g)
  n <- nrow(g)
  patt <- (g !=0 ) + (t(g) != 0)
  patt[patt>1] <- 1L
  R <- INLA::inla.as.sparse(
    Diagonal(n, 1L + rowSums(patt)) -(patt))
  ii <- INLA:::inla.qreordering(R)
  ##    return(ii)
  oR <- R[ii$ireordering, ii$ireordering]
  ##    print(oR)
  ij <- which(oR@i > oR@j)
  nnzq <- length(ij)
  if(is.null(theta))
    theta <- rnorm(n + nnzq)
  L <- Diagonal(n, exp(theta[1:n]))
  if(nnzq==0)
    return(crossprod(L))
  orL <- sparseMatrix(
    i = oR@i[ij] + 1L,
    j = oR@j[ij] + 1L,
    x = theta[n + 1:nnzq],
    dims = c(n, n)
  ) + L
  ##    print(orL)
  sL <- t(chol(as.matrix(oR)))
  lfi <- ((sL!=0) & (as.matrix(oR)==0)) * 1L
  if(sum(lfi)>0) {
##    cat("fill-in", sum(lfi), 'values\n')
    oL <- fiL(orL, lfi)
  }
  ##return(orL)
  return(tcrossprod(orL)[ii$reordering, ii$reordering])
}
