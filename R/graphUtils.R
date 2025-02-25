#' Set a graph where each node represent a variable
#' and the edges its conditional distribution. See details.
#' @param ... list of formula used to define the edges.
#' @details
#' The terms in the formula do represent the nodes.
#' The `~` is taken as link.
#' @export
#' @examples
#' g1 <- corgraph(x ~ y, y ~ v, v ~ z, z ~ x)
#' g1
#' summary(g1)
#' plot(g1)
#' precision(g1)
corgraph <- function(...) {
  fch <- as.character(match.call())[-1]
  if(length(fch)<1)
    stop("Please provide a graph argument!")
  ch <- lapply(fch, function(x)
    as.character(as.formula(x)))
  nodesL <- unique(unlist(
    lapply(ch, function(x) x[2])))
  ## right side check, and collect terms
  m <- length(ch)
  terms.r <- vector("list", m)
  for(i in 1:m) {
    x <- gsub(" ", "", ch[[i]][3])
    schi <- strsplit(x, "-", fixed = TRUE)[[1]]
    ##    print(schi)
    schi <- unlist(strsplit(schi, "+", fixed = TRUE))
    ##  print(schi)
    if(schi[1]=="") schi <- schi[-1]
    terms.r[[i]] <- schi
  }
  nodesR <- unique(unlist(terms.r))
  nodes <- unique(c(nodesL, nodesR))
  nNodes <- length(nodes)
  graph <- matrix(
    0, m, length(nodesR),
    dimnames = list(nodesL, nodesR))
  for(i in 1:m) {
    jj <- pmatch(terms.r[[i]], nodesR)
    graph[i, jj] <- 1
  }
  class(fch) <- 'corgraph'
  attr(fch, 'nodes') <- nodes
  attr(fch, 'graph') <- graph
  return(fch)
}
#' @export
print.corgraph <- function(x, ...) {
  cat("Model corgraph for",
      length(attr(x, 'nodes')), "variables",
      "using", sum(attr(x, 'graph')), "edges.\n")
}
#' @export
summary.corgraph <- function(object, ...) {
  attr(object, "graph")
}
#' @export
dim.corgraph <- function(x, ...) {
  c(nodes=length(attr(x, 'nodes')),
    edges=sum(attr(x, 'graph')))
}
setMethod(
  "edges",
  "corgraph",
  function(object, which, ...) {
    nodes <- attr(object, "nodes")
    stopifnot(!is.null(nodes))
    graph <- attr(object, "graph")
    stopifnot(!is.null(graph))
    m <- nrow(graph)
    edgl <- vector("list", m)
    er <- lapply(1:ncol(graph), function(i)
      colnames(graph)[graph[i, ]!=0])
    names(edgl) <- paste0(
      rownames(graph), "~",
      sapply(er, paste, collapse = "+"))
    for(i in 1:m) {
      edgl[[i]] <- list(
        n = sum(graph[i,]!=0),
        edges = er[[i]],
        weights = rep(1, length(er[[i]]))
      )
      edgl[[i]]$term <- edgl$edges
    }
    return(edgl)
  }
)
#' @export
setMethod(
  "plot",
  "corgraph",
  function(x, y, ...) {
    edgl <- edges(x)
    nodes <- names(edgl)
    gr <- graph::graphNEL(
      nodes = nodes,
      edgeL = edgl,
      edgemode='directed')
    print(gr)
    plot(gr, ...)
  }
)
#' The Laplacian of a matrix
#' @export
Laplacian.matrix <- function(graph) {
  if(inherits(graph, "matrix")) {
    A <- 1 - is.zero(graph)
    if(any(A!=t(A)))
      warning("Not symmetric!")
    L <- diag(rowSums(A)) - A
  } else {
    Laplacian.default(graph)
  }
}
#' The Laplacian method for 'corgraph'
#' @export
Laplacian.corgraph <- function(graph) {
  ne <- dim(graph)
  nodes <- attr(graph, "nodes")
  graph <- attr(graph, "graph")
  L <- matrix(
    0.0, ne[1], ne[1],
    dimnames = list(nodes, nodes)
    )
  for(i in 1:nrow(graph)) {
    ii <- pmatch(rownames(graph)[i], nodes)
    jj <- pmatch(colnames(graph)[graph[i, ]!=0], nodes)
    L[i, jj] <- (-1)
  }
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  return(L)
}
#' The precision method for 'corgraph'
#' @export
precision.corgraph <- function(x, ...) {
  Q <- Laplacian(x)
  n <- ncol(Q)
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
    ij <- which(lower.tri(Q) & !is.zero(Q))
    L <- diag(exp(theta[1:n]))
    L[ij] <- theta[-(1:n)]
    ll <- t(chol(Q + diag(1.0, n, n)))
    ifill <- which(is.zero(Q) & (!is.zero(ll)))
    if(length(ifill)>0) {
      L <- fiL(L, ifill)
    }
    Q <- tcrossprod(L)
  }
  return(Q)
}
#' The variance method for corgraph
#' @export
variance.corgraph <- function(x, ...) {
  Q <- precision(x, ...)
  return(chol2inv(chol(Q)))
}
#' Function to fill-in a Cholesky matrix
#' @param L the lower triangle of the Cholesky decomposition
#' @param lfi indicator of fill-in elements
#' @return a filled L matrix.
fiL <- function(L, lfi) {
  L <- as.matrix(L)
  ii <- row(L)[lfi]
  jj <- col(L)[lfi]
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
#' The is.zero.corgraph method
#' @export
is.zero.corgraph <- function(graph) {
  if(inherits(graph, "matrix"))
  Q <- precision(graph)
  is.zero.default(Q)
}
#' Evaluate the hessian of the KLD for graph model
#' around a base model.
#' @param graph model definition of a graphical model.
#' This can be either a matrix or a 'corgraph'.
#' @param base reference correlation matrix
#' or lower triangular parameters of a 'corgraph' model.
#' @param method the decomposition method used to
#' compute H^0.5 and H^(1/2).
#' @return list containing the hessian,
#' its 'square root', inverse 'square root' along
#' with the decomposition used
#' @examples
#' g <- corgraph(x ~ y+z, y ~ v, v ~ z)
#' ne <- dim(g)
#' gH0 <- graph2H(g, rep(-1, ne[2]))
#' ## alternatively
#' Q0 <- precision(g, theta = rep(c(0,-1), ne))
#' C0 <- cov2cor(solve(Q0))
#' all.equal(graph2H(g, C0), gH0)
#' @export
graph2H <- function(graph, base, method = c("eigen", "svd", "chol")) {
  method <- match.arg(method)
  graph2L <- function(theta) {
    L <- diag(n)
    L[lower.tri(L) & (!z0)] <- theta
    ifill <- which((!is.zero(l1)) & z0)
    if(length(ifill)>0) {
      L <- fiL(L, ifill)
    }
    return(L)
  }
  graph2C <- function(theta) {
    lQ <- graph2L(theta)
    V <- chol2inv(t(lQ))
    return(cov2cor(V))
  }
  Q0 <- Laplacian(graph)
  nEdges <- sum((!is.zero(Q0)) & lower.tri(Q0))
  z0 <- is.zero(Q0)
  n <- nrow(Q0)
  l1 <- t(chol(Q0 + diag(1.0, n, n)))
  if(inherits(base, "matrix")) {
    ## maybe optim() to get theta.base that
    ## give graph2C(theta.base) close to C0
    ## for now take L elements from C0
    C0 <- cov2cor(base)
    qq0 <- chol2inv(chol(C0))
    ll0 <- t(chol(qq0))
    for(i in 1:nrow(C0))
      ll0[i, ] <- ll0[i, ]/ll0[i, i]
    base <- ll0[lower.tri(ll0) & (!z0)]
  } else {
    stopifnot(length(base) == nEdges)
    C0 <- graph2C(base)
  }
  ## hessian uses graphpcor:::KLD10
  H <- hessian(function(x) KLD10(graph2C(x), C0),
               base)
  ## next bit follows mvtnorm:::rmvnorm()
  t0 <- sqrt(.Machine$double.eps)
  if(method == "eigen") {
    Hd <- eigen(H)
    if(!all(Hd$values >= (t0 * abs(Hd$values[1]))))
      warning("'H' is numerically not positive semidefinite")
    s <- sqrt(pmax(Hd$values, 0.0))
    h.5 <- t(Hd$vectors %*% (t(Hd$vectors) * s))
    hneg.5 <- t(Hd$vectors %*% (t(Hd$vectors) / s))
  }
  if(method == "svd") {
    Hd <- svd(H)
    if(any(Hd$d<(t0 * abs(Hd$d[1]))))
      warning("'H' is numerically not positive semidefinite")
    s <- sqrt(pmax(Hd$d, 0.0))
    h.5 <- t(Hd$v %*% (t(Hd$u) * s))
    hneg.5 <- t(Hd$v %*% (t(Hd$u) / s))
  }
  if(method == "chol") {
    Hd <- chol(H, pivot = TRUE)
    h.5 <- matrix(Hd[, order(attr(Hd, "pivot")), ], nrow(H))
    hn <- chol2inv(chol(H))
    hn.5 <- chol(Hn, pivot = TRUE)
    hneg.5 <- matrix(hn.5[, order(attr(hn.5, "pivot"))], nrow(H))
  }
  stopifnot(all.equal(H, tcrossprod(h.5)))
  list(H = H,
       h.5 = h.5,
       hneg.5 = hneg.5,
       Hdecomposition = Hd)
}
