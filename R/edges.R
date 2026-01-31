#' Retrieve edges of an object
#' @rdname edges
#' @export
edges <- function(object) {
  UseMethod("edges")
}
#' @describeIn edges
#' Default method for edges
#' @export
edges.default <- function(object) {
  stop("No 'edges' method for this class!")
}
#' @describeIn edges
#' Extract the edges of a `graphpcor`.
#' @param object graphpcor object
#' @export
edges.graphpcor <- function(object) {
  ne <- dim(object)
  nodes <- attr(object, "nodes")
  stopifnot(!is.null(nodes))
  stopifnot(ne[1]==length(nodes))
  Q1 <- Laplacian(object)
  edgl <- vector("list", ne[1])
  for(i in 1:ne[1]) {
    jj <- setdiff(which(!is.zero(Q1[i, ])), i)
    ni <- length(jj)
    if(ni>0) {
      edgl[[i]] <- list(
        n = ni,
        edges = nodes[jj],
        weights = rep(1.0, ni))
      edgl[[i]]$term <- jj
    }
  }
  names(edgl) <- nodes
  return(edgl)
}
#' @describeIn edges
#' Extract the edges of a `treepcor`
#' @importFrom methods new
#' @export
edges.treepcor <- function(object) {
  trm <- attr(object, "relationship")
  m <- ncol(trm)
  n <- nrow(trm)-m+1
  if(m>1) {
    trm <- trm[c(n+1:(m-1), 1:n), ]
    stopifnot(all(substr(rownames(trm)[1:(m-1)],1,1) == "p"))
  }
  stopifnot(all(substr(rownames(trm)[m:nrow(trm)],1,1) == "c"))
  stopifnot(all(substr(rownames(trm),1,1) %in% c("p", "c")))
  edgl <- vector("list", m + n)
  names(edgl) <- c(paste0("p", 1:m),
                   paste0("c", 1:n))
  for(i in 1:m) {
    w1 <- trm[, i] != 0
    edgl[[i]] <- list(
      n = sum(w1),
      edges = rownames(trm)[w1],
      weights = new("numeric", trm[w1, i])
    )
    edgl[[i]]$term <- edgl$edges
    edgl[[i]]$parent <- substr(edgl[[i]]$edges, 1, 1) == "p"
    edgl[[i]]$id <- new("integer", substring(edgl[[i]]$edges, 2))
    edgl[[i]]$sign <- ifelse(edgl[[i]]$weights<0, -1, 1)
  }
  return(edgl)
}
