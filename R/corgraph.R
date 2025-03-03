#' @describeIn corgraph
#' The `corgraph` is a corgraph where each node represents
#' a variable and each edge indicates a conditional distribution.
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
corgraph.formula <- function(...) {
  fch <- as.character(match.call())[-1]
  if(length(fch)<1)
    stop("Please provide an argument!")
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
  grel <- matrix(
    0, m, length(nodesR),
    dimnames = list(nodesL, nodesR))
  for(i in 1:m) {
    jj <- pmatch(terms.r[[i]], nodesR)
    grel[i, jj] <- 1
  }
  class(fch) <- 'corgraph'
  attr(fch, 'nodes') <- nodes
  attr(fch, 'relationship') <- grel
  return(fch)
}
#' @describeIn corgraph
#' Build a `corgraph` from a matrix
#' @export
corgraph.matrix <- function(x) {
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
  return(do.call(what = 'corgraph',
                 args = lapply(argl, as.formula)))
}
#' @describeIn corgraph
#' The print method for `corgraph`
#' @export
print.corgraph <- function(x, ...) {
  cat("A corgraph for",
      length(attr(x, 'nodes')), "variables",
      "using", sum(attr(x, 'graph')), "edges.\n")
}
#' @describeIn corgraph
#' The summary method for `corgraph`
#' @export
summary.corgraph <- function(object, ...) {
  attr(object, "relationship")
}
#' @describeIn corgraph
#' The dim method for `corgraph`
#' @export
dim.corgraph <- function(x, ...) {
  c(nodes=length(attr(x, 'nodes')),
    edges=sum(attr(x, 'relationship')))
}
#' @describeIn corgraph
#' The plot method for `corgraph`
#' @export
setMethod(
  "plot",
  "corgraph",
  function(x, y, ...) {
    ne <- dim(x)
    nodes <- attr(x, "nodes")
    stopifnot(!is.null(nodes))
    stopifnot(ne[1]==length(nodes))
    L <- Laplacian(x)
    edgl <- vector("list", ne[1])
    for(i in 1:ne[1]) {
      jj <- setdiff(which(!is.zero(L[i, ])), i)
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
    gr <- graph::graphNEL(
      nodes = nodes,
      edgeL = edgl,
      edgemode = 'undirected')
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
#' @describeIn corgraph
#' The Laplacian method for `corgraph`
#' @export
Laplacian.corgraph <- function(graph) {
  ne <- dim(graph)
  nodes <- attr(graph, "nodes")
  grel <- attr(graph, "relationship")
  L <- matrix(
    0.0, ne[1], ne[1],
    dimnames = list(nodes, nodes)
  )
  for(i in 1:nrow(grel)) {
    ii <- pmatch(rownames(grel)[i], nodes)
    jj <- pmatch(colnames(grel)[grel[i, ]!=0], nodes)
    L[i, jj] <- (-1)
  }
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  return(L)
}
#' @describeIn corgraph
#' Build the unite diagonal lower triangle matrix
setMethod(
  "chol",
  "corgraph",
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
#' @describeIn corgraph
#' The variance method for 'corgraph'
#' @export
variance.corgraph <- function(x, ...) {
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
  L <- getMethod('chol', 'corgraph')(
    x, theta = theta[-(1:ne[1])])
  V <- chol2inv(L)
  si <- exp(theta[1:ne[1]]) / sqrt(diag(V))
  V <- diag(si) %*% V %*% diag(si)
  return(V)
}
#' @describeIn corgraph
#' The precision method for 'corgraph'
#' @export
precision.corgraph <- function(x, ...) {
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
#' Evaluate the hessian of the KLD for a `corgraph`
#' correlation model around a base model.
#' @param corgraph model definition of a graphical model.
#' This can be either a matrix or a 'corgraph'.
#' @param base either a reference correlation matrix
#' or as a parameter reference for a 'corgraph' model.
#' @param decomposition character to specify which
#' decomposition method to use to compute H^0.5 and H^(1/2).
#' @param ... additional arguments passed to [numDeriv::hessian()]
#' @return list containing the hessian,
#' its 'square root', inverse 'square root' along
#' with the decomposition used
#' @examples
#' g <- corgraph(x1 ~ x2 + x3, x2 ~ x4, x3 ~ x4)
#' ne <- dim(g)
#' gH0 <- hessian(g, rep(-1, ne[2]))
#' ## alternatively
#' Q0 <- precision(g, theta = rep(c(0,-1), ne))
#' C0 <- cov2cor(solve(Q0))
#' all.equal(hessian(g, C0), gH0)
#' @export
hessian.corgraph <- function(corgraph, base, decomposition = c("eigen", "svd", "chol"), ...) {
  decomposition <- match.arg(decomposition)
  Q0 <- Laplacian(corgraph)
  nEdges <- sum((!is.zero(Q0)) & lower.tri(Q0))
  z0 <- is.zero(Q0)
  n <- nrow(Q0)
  l1 <- t(chol(Q0 + diag(1.0, n, n)))
  if(inherits(base, "matrix")) {
    ## maybe optim() to get theta.base that
    ## give it close to C0?
    ## For now check the elements of L from C0^{-1}
    C0 <- cov2cor(base)
    qq0 <- chol2inv(chol(C0))
    ll0 <- t(chol(qq0))
    c0.ok <- all(which(abs(ll0)>sqrt(.Machine$double.eps)) %in%
                   which(abs(l1)>0))
    if(!c0.ok) {
      stop("Provided base correlation not in the corgraph model class!")
    }
    for(i in 1:nrow(C0))
      ll0[i, ] <- ll0[i, ]/ll0[i, i]
    base <- ll0[lower.tri(ll0) & (!z0)]
  } else {
    stopifnot(length(base) == nEdges)
    C0 <- variance(corgraph, theta =base)
  }
  ## hessian uses graphpcor:::KLD10
  H <- numDeriv::hessian(
    function(x) KLD10(variance(corgraph, theta = x), C0),
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
#' The `cgeneric` method for `corgraph` uses [cgeneric_corgraph()]
#' @export
cgeneric.corgraph <- function(...) {
  args <- list(...)
  do.call(what = 'cgeneric_corgraph',
          args = args)
}
