#' Function to check a graph
#' @param graph the graph
graph_check <- function(graph) {
  nlinks <- length(graph)
  gchar <- gsub(" ", "", as.character(graph))
  stopifnot(length(unique(grep("~", gchar))) == nlinks)
  chs <- sapply(gchar, function(x)
    strsplit(x, split = "~")[[1]])
  r <- all(substr(chs[1, ], 1, 1) == "c") &
    all(substr(chs[2, ], 1, 1) == "c")
  ij <- as.integer(substring(chs, 2))
  r <- r&all(!is.na(ij))
  attr(r, 'chs') <- chs
  attr(r, 'ij') <- ij
  return(r)
}
#' Function to create link indexes from a graph
graph_elements <- function(graph) {
  test <- graph_check(graph)
  stopifnot(test)
  ii <- as.integer(substring(attr(test, "chs")[1, ], 2))
  jj <- as.integer(substring(attr(test, "chs")[2, ], 2))
  if(all(ii<jj)) {
    ret <- list(ii = jj, jj = ii)
  } else {
    stopifnot(all(ii>jj))
    ret <- list(ii = ii, jj = jj)
  }
  return(ret)
}
#' Precision structure (as discrete Laplacian)
#' @param ij output of graph_elements
#' @export
graph_Laplacian <- function(graph) {
  ij <- graph_elements(graph)
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
  attr(q, "ii") <- ij$ii
  attr(q, "jj") <- ij$jj
  return(q)
}
#' Precision and Cholesky fill-in indexes
graph_qchol_index <- function(graph) {
  Lap <- graph_Laplacian(graph)
  n <- nrow(Lap)
  ret <- list(
    Q = Lap,
    n = n,
    ii = attr(Lap, "ii"),
    jj = attr(Lap, "jj")
    )
  qnz <- Lap!=0
  ret$iq <- which(qnz)
  ret$ilq <- which(
    qnz & lower.tri(q, diag = TRUE))
  ll <- t(chol(Lap + diag(n)))
  ret$ifil <- intersect(ret$ilq, which(ll!=0))
  return(ret)
}
