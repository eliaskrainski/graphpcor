#' @describeIn graph
#' The `graph` is a graph where each node represents
#' a variable and each edge indicates a conditional distribution.
#' @param ... list of formula used to define the edges.
#' @details
#' The terms in the formula do represent the nodes.
#' The `~` is taken as link.
#' @export
#' @examples
#' g1 <- graph(x ~ y, y ~ v, v ~ z, z ~ x)
#' g1
#' summary(g1)
#' plot(g1)
#' precision(g1)
graph.formula <- function(...) {
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
  class(fch) <- 'graph'
  attr(fch, 'nodes') <- nodes
  attr(fch, 'graph') <- graph
  return(fch)
}
#' @describeIn graph
#' Build a graph from a matrix
#' @export
graph.matrix <- function(x) {
  #  if(inherits(list(...)[[1]], "matrix")) {
  #    x <- list(...)[[1]]
  stopifnot(all.equal(x, t(x)))
  ne <- c(nrow(x), NA)
  iz <- is.zero(x)
  ne[2] <- (sum(!iz)-ne[1])/2
  vnams <- rownames(x)
  if(is.null(vnams)) {
    vnams <- letters[1:ne[1]]
  }
  argl <- lapply(1:(ne[1]-1), function(i) {
    jj <- intersect((i+1):ne[1], which(!iz[i, ]))
    paste(vnams[i], "~",
          paste(vnams[jj], collapse = " + "))
  })
  return(do.call(what = 'graph',
                 args = lapply(argl, as.formula)))
}
#' @describeIn graph
#' The print method for `graph`
#' @export
print.graph <- function(x, ...) {
  cat("Model graph for",
      length(attr(x, 'nodes')), "variables",
      "using", sum(attr(x, 'graph')), "edges.\n")
}
#' @describeIn graph
#' The summary method for `graph`
#' @export
summary.graph <- function(object, ...) {
  attr(object, "graph")
}
#' @describeIn graph
#' The dim method for `graph`
#' @export
dim.graph <- function(x, ...) {
  c(nodes=length(attr(x, 'nodes')),
    edges=sum(attr(x, 'graph')))
}
#' @describeIn graph
#' The edges method for `graph`
setMethod(
  "edges",
  "graph",
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
#' @describeIn graph
#' The plot method for `graph`
#' @export
setMethod(
  "plot",
  "graph",
  function(x, y, ...) {
    edgl <- edges(x)
    nodes <- names(edgl)
    gr <- graph::graphNEL(
      nodes = nodes,
      edgeL = edgl,
      edgemode='directed')
    plot(gr, ...)
  }
)
#' @describeIn Laplacian
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
#' @describeIn graph
#' The Laplacian method for `graph`
#' @export
Laplacian.graph <- function(graph) {
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
#' @describeIn graph
#' Build the unite diagonal lower triangle matrix
setMethod(
  "chol",
  "graph",
  function(x, ...) {
    ne <- dim(x)
    G <- Laplacian(x)
    idx <- which(lower.tri(G) & (!is.zero(G)))
    L <- diag(ne[1])
    stopifnot(length(idx)==ne[2])
    args <- list(...)
    stopifnot(length(args$theta)==ne[2])
    L[idx] <- args$theta
    ll <- t(chol(G + diag(ne[1])))
    ifill <- which(is.zero(G) & (!is.zero(ll)))
    if(length(ifill)>0) {
      L <- fiL(L, ifill)
    }
    return(t(L))
  }
)
#' @describeIn graph
#' The variance method for 'graph'
#' @export
variance.graph <- function(x, ...) {
  mc <- list(...)
  nargs <- names(mc)
  if(!any(nargs == "theta")) {
    stop("Please provide 'theta'!")
  }
  theta <- mc$theta
  ne <- dim(x)
  Q <- Laplacian(x)
  stopifnot(ne[1]==nrow(Q))
  stopifnot((2*ne[2])==(sum(!is.zero(Q))-ne[1]))
  if(length(theta)==ne[2]) {
    theta <- c(rep(0.0, ne[1]), theta)
  } else {
    stopifnot(length(theta)==sum(ne))
  }
  L <- getMethod('chol', 'graph')(
    x, theta = theta[-(1:ne[1])])
  V <- chol2inv(L)
  si <- exp(theta[1:ne[1]]) / sqrt(diag(V))
  V <- diag(si) %*% V %*% diag(si)
  return(V)
}
#' @describeIn graph
#' The precision method for 'graph'
#' @export
precision.graph <- function(x, ...) {
  ne <- dim(x)
  Q <- Laplacian(x)
  stopifnot(ne[1]==nrow(Q))
  stopifnot((2*ne[2])==(sum(!is.zero(Q))-ne[1]))
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    return(chol2inv(chol(
      variance(x, ...))
    ))
  } else {
    warning("missing `theta`, returning Laplacian!")
  }
  return(Q)
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
#' @describeIn is.zero
#' The is.zero.graph method
#' @export
is.zero.graph <- function(graph) {
  if(inherits(graph, "matrix"))
    Q <- precision(graph)
  is.zero.default(Q)
}
#' Evaluate the hessian of the KLD for graph model
#' around a base model.
#' @param graph model definition of a graphical model.
#' This can be either a matrix or a 'graph'.
#' @param base either a reference correlation matrix
#' or as a parameter reference for as 'graph' model.
#' @param decomposition character to specify which
#' decomposition method to use to compute H^0.5 and H^(1/2).
#' @param ... additional arguments passed to [numDeriv::hessian()]
#' @return list containing the hessian,
#' its 'square root', inverse 'square root' along
#' with the decomposition used
#' @examples
#' g <- graph(x ~ y+z, y ~ v, v ~ z)
#' ne <- dim(g)
#' gH0 <- hessian(g, rep(-1, ne[2]))
#' ## alternatively
#' Q0 <- precision(g, theta = rep(c(0,-1), ne))
#' C0 <- cov2cor(solve(Q0))
#' all.equal(hessian(g, C0), gH0)
#' @export
hessian.graph <- function(graph, base, decomposition = c("eigen", "svd", "chol"), ...) {
  decomposition <- match.arg(decomposition)
  Q0 <- Laplacian(graph)
  nEdges <- sum((!is.zero(Q0)) & lower.tri(Q0))
  z0 <- is.zero(Q0)
  n <- nrow(Q0)
  l1 <- t(chol(Q0 + diag(1.0, n, n)))
  if(inherits(base, "matrix")) {
    ## maybe optim() to get theta.base that
    ## give graph2C(theta.base) close to C0?
    ## For now check the elements of L from C0^{-1}
    C0 <- cov2cor(base)
    qq0 <- chol2inv(chol(C0))
    ll0 <- t(chol(qq0))
    c0.ok <- all(which(abs(ll0)>sqrt(.Machine$double.eps)) %in%
                   which(abs(l1)>0))
    if(!c0.ok) {
      stop("Provided base correlation not in the graph model class!")
    }
    for(i in 1:nrow(C0))
      ll0[i, ] <- ll0[i, ]/ll0[i, i]
    base <- ll0[lower.tri(ll0) & (!z0)]
  } else {
    stopifnot(length(base) == nEdges)
    C0 <- variance(graph, theta =base)
  }
  ## hessian uses graphpcor:::KLD10
  H <- numDeriv::hessian(
    function(x) KLD10(variance(graph, theta = x), C0),
    base, ...)
  ## next bit follows mvtnorm:::rmvnorm()
  t0 <- sqrt(.Machine$double.eps)
  if(decomposition == "eigen") {
    Hd <- eigen(H)
    if(!all(Hd$values >= (t0 * abs(Hd$values[1]))))
      warning("'H' is numerically not positive semidefinite")
    s <- sqrt(pmax(Hd$values, 0.0))
    h.5 <- t(Hd$vectors %*% (t(Hd$vectors) * s))
    hneg.5 <- t(Hd$vectors %*% (t(Hd$vectors) / s))
  }
  if(decomposition == "svd") {
    Hd <- svd(H)
    if(any(Hd$d<(t0 * abs(Hd$d[1]))))
      warning("'H' is numerically not positive semidefinite")
    s <- sqrt(pmax(Hd$d, 0.0))
    h.5 <- t(Hd$v %*% (t(Hd$u) * s))
    hneg.5 <- t(Hd$v %*% (t(Hd$u) / s))
  }
  if(decomposition == "chol") {
    Hd <- chol(H, pivot = TRUE)
    h.5 <- matrix(Hd[, order(attr(Hd, "pivot")), ], nrow(H))
    hn <- chol2inv(chol(H))
    hn.5 <- chol(Hn, pivot = TRUE)
    hneg.5 <- matrix(hn.5[, order(attr(hn.5, "pivot"))], nrow(H))
  }
  stopifnot(all.equal(H, tcrossprod(h.5)))
  attr(H, "base") <- base
  attr(H, "h.5") < h.5
  attr(H, "hneg.5") <- hneg.5
  attr(H, "decomposition") <- Hd
  return(H)
}
#' @describeIn cgeneric
#' The `cgeneric` method for `graph` uses [cgeneric_pcgraph()]
#' @export
cgeneric.graph <- function(...) {
  args <- list(...)
  args$graph <- args$model
  args$model <- NULL
  do.call(what = 'cgeneric_pcgraph',
          args = args)
}
