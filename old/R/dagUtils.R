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
##  stopifnot(n == max(ii,jj))
  nnzq <- length(ii)
  if(is.null(theta))
    theta <- rnorm(n + nnzq)
  L <- diag(exp(theta[1:n]))
  if(nnzq==0)
    return(INLA::inla.as.sparse(L))
  retL <- as.matrix(sparseMatrix(
    i = ii,
    j = jj,
    x = theta[n + 1:nnzq],
    dims = c(n, n)
  )) + L
  sL <- t(chol(R + diag(n)))
  lfi <- which((sL!=0) & (R==0))
  if(length(lfi)>0) {
    retL <- fiL(retL, lfi)
  }
  return(INLA::inla.as.sparse(retL))
}
#' Function to check a graph
#' @param graph the graph
graph_check <- function(graph) {
  nlinks <- sum(!sapply(graph, is.null))
  if(nlinks>0) {
    gchar <- gsub(" ", "", as.character(graph))
    stopifnot(length(unique(grep("~", gchar))) == nlinks)
    chs <- sapply(gchar, function(x)
      strsplit(x, split = "~")[[1]])
    r <- all(substr(chs[1, ], 1, 1) == "c") &
      all(substr(chs[2, ], 1, 1) == "c")
    ij <- as.integer(substring(chs, 2))
    #    print(ij)
    r <- r&all(!is.na(ij))
    n <- max(ij)
  } else {
    r <- TRUE
    chs <- matrix("", 2, 0)
    ij <- integer(0)
  }
  attr(r, 'n') <- n
  attr(r, 'chs') <- chs
  attr(r, 'ij') <- ij
  return(r)
}
#' Function to create link indexes from a graph
graph_elements <- function(graph) {
  test <- graph_check(graph)
  stopifnot(test)
  ij <- matrix(attr(test, "ij"), nrow = 2)
  ii <- ij[1, ]
  jj <- ij[2, ]
  if(all(ii<jj)) {
    ret <- list(ii = jj, jj = ii)
  } else {
    stopifnot(all(ii>jj))
    ret <- list(ii = ii, jj = jj)
  }
  attr(ret, "n") <- attr(test, "n")
  return(ret)
}
#' Precision structure (as discrete Laplacian)
#' @param ij output of graph_elements
#' @export
graph_Laplacian <- function(graph) {
  ij <- graph_elements(graph)
  n <- attr(ij, "n") ##max(ij$ii, ij$jj)
  q <- matrix(0, n, n)
  if(length(ij$ii)>0) {
    for(k in 1:length(ij$ii)) {
      i <- ij$ii[k]
      j <- ij$jj[k]
      q[i, i] <- q[i, i] +1
      q[j, j] <- q[j, j] +1
      q[i, j] <- q[i, j] -1
      q[j, i] <- q[j, i] -1
    }
  }
  attr(q, "ii") <- ij$ii
  attr(q, "jj") <- ij$jj
  return(q)
}
#' Precision and Cholesky fill-in indexes
graph_qchol_index <- function(graph) {
  Lap <- graph_Laplacian(graph_Laplacian(graph))
  n <- nrow(Lap)
  Lap <- Lap + diag(n)
  ret <- list(
    Q = Lap,
    n = n,
    ii = attr(Lap, "ii"),
    jj = attr(Lap, "jj")
  )
  qnz <- Lap!=0
  ret$iq <- which(qnz)
  ret$ilq <- which(
    qnz & lower.tri(Lap, diag = TRUE))
  ret$iuq <- which(
    qnz & upper.tri(Lap, diag = TRUE))
  ll <- t(chol(Lap + diag(n)))
  ret$ifil <- setdiff(which(ll!=0), ret$ilq)
  return(ret)
}
