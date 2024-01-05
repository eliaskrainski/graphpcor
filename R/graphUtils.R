#' Function to check a graph
#' @param graph the graph
graph_check <- function(graph) {
  nlinks <- length(graph)
  gchar <- gsub(" ", "", as.character(graph))
  stopifnot(length(unique(grep("~", gchar)))==nlinks)
  gsch <- sapply(gchar, function(x)
    strsplit(x, split = "~")[[1]])
  stopifnot(all(substr(gsch[1, ], 1, 1) == "c"))
  stopifnot(all(substr(gsch[2, ], 1, 1) == "c"))
  r <- TRUE
  attr(r, "gsch") <- gsch
  return(r)
}
#' Function to create link indexes from a graph
graph_elements <- function(graph) {
  test <- graph_check(graph)
  stopifnot(test)
  ii <- as.integer(substring(attr(test, "gsch")[1, ], 2))
  jj <- as.integer(substring(attr(test, "gsch")[2, ], 2))
  return(list(ii = ii, jj = jj))
}
#' Precision structure (as discrete Laplacian)
#' @param ij output of graph_elements
graph_Laplacian <- function(ij) {
  n <- max(ij$ii, ij$jj)
  q <- matrix(0, n, n)
  for(k in 1:length(ij$ii)) {
    i <- ij$ii[k]
    j <- ij$jj[k]
    q[i, i] <- q[i, i] +1
    q[j, j] <- q[j, j] +1
    q[i, j] <- q[i, j] -1
    q[j, i] <- q[j, i] -1
  }
  return(q)
}
#' Precision and Cholesky fill-in indexes
graph_qchol_index <- function(graph) {
  ij <- graph_elements(graph)
  n <- max(ij$ii, ij$jj)
  q <- graph_Laplacian(ij) + diag(n)
  ret <- list(n = n, i = ij$ii, j = ij$jj)
  qnz <- q!=0
  ret$qii <- which(qnz)
  ret$qiiu <- which(qnz & upper.tri(q, diag = TRUE))
  l <- chol(q)
  ll <- abs(l*1000)
  ret$lii <- which(ll!=0)
  return(ret)
}
