#' Make a graph from a matrix
#' @param G the link matrix
matrix2graph <- function(G) {
  iG <- (G!=0) & upper.tri(G)
  i <- row(G)[iG]
  j <- col(G)[iG]
  nlinks <- length(i)
  if(nlinks<1) return(NULL)
  r <- lapply(1:nlinks, function(k)
    as.formula(paste0("c", i[k], "~c", j[k])))
  return(r)
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
  Lap <- graph_Laplacian(graph)
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
