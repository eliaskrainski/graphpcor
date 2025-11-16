#' @describeIn graphpcor
#' A `graphpcor` is a graph where a node represents
#' a variable and an edge a conditional distribution.
#' @param ... list of formula used to define the edges.
#' @details
#' The terms in the formula do represent the nodes.
#' The `~` is taken as link.
#' @importFrom stats as.formula
#' @export
#' @examples
#' g1 <- graphpcor(x ~ y, y ~ v, v ~ z, z ~ x)
#' g1
#' summary(g1)
#' plot(g1)
#' prec(g1)
graphpcor.formula <- function(...) {
  fch <- as.character(match.call())[-1]
  m <- length(fch)
  if(m<1)
    stop("Please provide an argument!")
  ch <- lapply(fch, function(x)
    as.character(as.formula(x)))
  ## right side check, and collect terms
  terms.r <- lapply(ch, function(x) {
    x <- gsub(" ", "", x[3])
    schi <- strsplit(x, "-", fixed = TRUE)[[1]]
    schi <- unlist(strsplit(schi, "+", fixed = TRUE))
    if(schi[1]=="") schi <- schi[-1]
    return(schi)
  })
  nodesL <- sapply(ch, function(x) x[2])
  allNodes <- unique(unlist(lapply(1:m, function(i)
    c(nodesL[i], terms.r[[i]]))))
  nNodes <- length(allNodes)
  ## graph
  grel <- matrix(0, nNodes, nNodes,
                 dimnames = list(allNodes, allNodes))
  for(i in 1:m) {
    ii <- pmatch(nodesL[i], allNodes)
    jj <- pmatch(terms.r[[i]], allNodes)
    if((length(ii)>0) & (length(jj)>0)) {
      grel[ii, jj] <- 1
      grel[jj, ii] <- 1
    }
  }
  class(fch) <- 'graphpcor'
  attr(fch, 'nodes') <- allNodes
  attr(fch, 'graph') <- grel
  return(fch)
}
#' @describeIn graphpcor
#' Build a `graphpcor` from a matrix
#' @importFrom stats as.formula
#' @importFrom INLAtools is.zero
#' @export
graphpcor.matrix <- function(...) {
  x <- list(...)[[1]]
  stopifnot(all.equal(x, t(x)))
  ne <- c(nrow(x), NA)
  iz <- is.zero(x, tol = 1e-9)
  ne[2] <- sum(lower.tri(iz) & (!iz))
  vnams <- rownames(x)
  if(is.null(vnams)) {
    vnams <- letters[1:ne[1]]
  }
  argl <- list()
  adde <- 0
  for(i in 1:(ne[1]-1)) {
    jj <- intersect(
      (i+1):ne[1],
      which(!iz[i, ]))
    if(length(jj)>0) {
      argl[[i]] <- paste(
        vnams[i], "~",
        paste(vnams[jj], collapse = " + "))
      adde <- adde + length(jj)
    }
  }
  stopifnot(adde==ne[2])
  return(do.call(what = 'graphpcor',
                 args = lapply(argl, as.formula)))
}
#' @describeIn graphpcor
#' The print method for `graphpcor`
#' @param x graphpcor
#' @export
print.graphpcor <- function(x, ...) {
  n <- length(attr(x, 'nodes'))
  g <- !is.zero(attr(x, 'graph'))
  cat("A graphpcor for", n, "variables",
      "with", sum(lower.tri(g) & g), "edges.\n")
}
#' @describeIn graphpcor
#' The summary method for `graphpcor`
#' @param object graphpcor
#' @export
summary.graphpcor <- function(object, ...) {
  attr(object, "graph")
}
#' @describeIn graphpcor
#' The dim method for `graphpcor`
#' @param x graphpcor
#' @export
dim.graphpcor <- function(x, ...) {
  c(nodes=length(attr(x, 'nodes')),
    edges=sum(attr(x, 'graph'))/2)
}
#' @describeIn graphpcor
#' Extract the edges of a `graphcor` to be used for plot
#' @param object graphpcor object
#' @param which not used
#' @importFrom graph edges
#' @export
setMethod(
  "edges",
  "graphpcor",
  function(object, which, ...) {
    ne <- dim(object)
    nodes <- attr(object, "nodes")
    stopifnot(!is.null(nodes))
    stopifnot(ne[1]==length(nodes))
    L <- Laplacian(object)
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
    return(edgl)
  }
)
#' @describeIn graphpcor
#' The plot method for `graphpcor`
#' @param x graphpcor
#' @param y graphpcor
#' @importFrom methods getMethod
#' @export
setMethod(
  "plot",
  "graphpcor",
  function(x, y, ...) {
    ne <- dim(x)
    nodes <- attr(x, "nodes")
    stopifnot(!is.null(nodes))
    stopifnot(ne[1]==length(nodes))
    edgl <- edges(x, which)
    gr <- graph::graphNEL(
      nodes = nodes,
      edgeL = edgl,
      edgemode = 'undirected')

    mc <- lapply(
      match.call(expand.dots = TRUE)[-1],
      eval)
    nargs <- names(mc)
    nattr <- list(color = {
      if(any(nargs=="color")) mc$color
      else rep("blue", ne[1])
    },
    fillcolor = {
      if(any(nargs == "fillcolor"))
        mc$fillcolor
      else rep("lightblue", ne[1])
    },
    shape = {
      if(any(nargs == "shape"))
        mc$shape
      else
        rep("circle", ne[1])
    },
    height = {
      if(any(nargs == "height"))
        mc$height
      else
        rep(0.5, ne[1])
    },
    width = {
      if(any(nargs == "width"))
        mc$width
      else
        rep(1.5, ne[1])
    },
    fontsize = {
      if(any(nargs == "fontsize"))
        mc$fontsize
      else
        rep(14, ne[1])
    }
    )
    for(i in 1:length(nattr))
      names(nattr[[i]]) <- nodes

    ag <- Rgraphviz::agopen(gr, "", nodeAttrs = nattr)

    for(k in 1:length(ag@AgEdge)) {
      ag@AgEdge[[k]]@color <- "red"
    }
    getMethod("plot", "Ragraph")(ag)
  }
)
#' @describeIn Laplacian
#' The Laplacian of a matrix
#' @param graph an object that inherits a matrix class
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
#' @describeIn graphpcor
#' The Laplacian method for a `graphpcor`
#' @param graph graphpcor object, see [`graphpcor`].
#' @export
Laplacian.graphpcor <- function(graph) {
  ne <- dim(graph)
  nodes <- attr(graph, "nodes")
  L <- -attr(graph, "graph")
  diag(L) <- -colSums(L)
  return(L)
}
#' @describeIn graphpcor
#' The `vcov` method for a `graphpcor`
#' @importFrom methods getMethod
#' @export
setMethod(
  "vcov",
  "graphpcor",
  function(object, ...) {
    mc <- list(...)
    nargs <- names(mc)
    if(!any(nargs == "theta")) {
      stop("Please provide 'theta'!")
    }
    theta <- mc$theta
    ne <- dim(object)
    G <- Laplacian(object)
    stopifnot(ne[1]==nrow(G))
    stopifnot((2*ne[2])==(sum(!is.zero(G))-ne[1]))
    if(length(theta)==ne[2]) {
      theta <- c(rep(0.0, ne[1]), theta)
    } else {
      stopifnot(length(theta)==sum(ne))
    }
    itheta <- which(lower.tri(G) & (!is.zero(G)))
    L <- t(theta2L(
      theta = theta[-(1:ne[1])],
      p = ne[1],
      parametrization = "itp",
      itheta = itheta,
      d0 = ne[1]:1))
    V <- chol2inv(L)
    si <- exp(theta[1:ne[1]]) / sqrt(diag(V))
    V <- diag(si) %*% V %*% diag(si)
    return(V)
  }
)
#' @describeIn graphpcor
#' The precision method for 'graphpcor'
#' @param model graphpcor model object
#' @importFrom Matrix Matrix forceSymmetric
#' @importFrom INLAtools Sparse prec
#' @export
prec.graphpcor <- function(model, ...) {
  ne <- dim(model)
  Q <- Laplacian(model)
  stopifnot(ne[1]==nrow(Q))
  stopifnot((2*ne[2])==(sum(!is.zero(Q))-ne[1]))
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
    if(length(theta)==ne[2]) {
      theta <- c(rep(0.0, ne[1]), theta)
    } else {
      stopifnot(length(theta)==sum(ne))
    }
    V <- vcov(model, theta = theta)
    Q <- chol2inv(chol(V))
    Q[is.zero(Q)] <- 0
  } else {
    warning("missing `theta`, returning Laplacian!")
  }
  return(forceSymmetric(
    Sparse(Matrix(Q),
           zeros.rm = TRUE)))
}
#' Evaluate the hessian of the KLD for a `graphpcor`
#' correlation model around a base model.
#' @param func model definition of a graphical model.
#' This can be either a matrix or a 'graphpcor'.
#' @param x either a reference correlation matrix
#' or a numeric vector with the parameters for the
#' reference 'graphpcor' model.
#' @param method see [numDeriv::hessian()]
#' @param method.args see [numDeriv::hessian()]
#' @param ... use to pass the decomposition method,
#' as a character to specify which one is to be used
#' to compute H^0.5 and H^(1/2).
#' @return list containing the hessian,
#' its 'square root', inverse 'square root' along
#' with the decomposition used
#' @examples
#' g <- graphpcor(x1 ~ x2 + x3, x2 ~ x4, x3 ~ x4)
#' ne <- dim(g)
#' gH0 <- hessian(g, rep(-1, ne[2]))
#' ## alternatively
#' C0 <- vcov(g, theta = rep(c(0,-1), ne))
#' all.equal(hessian(g, C0), gH0)
#' @importFrom stats cov2cor
#' @importFrom numDeriv hessian
#' @importFrom INLAtools is.zero
#' @export
hessian.graphpcor <- function(
    func,
    x,
    method = "Richardson",
    method.args = list(),
    ...) {
    decomposition <- c("svd", "eigen", "chol")
    if(is.null(list(...)$decomposition)) {
      decomposition <- "svd"
    } else {
      decomposition <- match.arg(
        list(...)$decomposition,
        decomposition)
    }
    Q0 <- Laplacian(func)
    nEdges <- sum((!is.zero(Q0, tol = 1e-9)) & lower.tri(Q0))
    z0 <- is.zero(Q0, tol = 1e-9)
    n <- nrow(Q0)
    l1 <- t(chol(Q0 + diag(1.0, n, n)))
    if(inherits(x, "Matrix")) {
      x <- as.matrix(x)
    }
    if(inherits(x, "matrix")) {
      ## maybe find theta that gives it close to C0?
      ## For now check the elements of L from C0^{-1}
      C0 <- cov2cor(x)
      qq0 <- chol2inv(chol(C0))
      ll0 <- t(chol(qq0))
      for(i in 1:n)
        ll0[i, ] <- (n+1-i)*ll0[i, ]/ll0[i, i]
      c0.ok <- all(which(abs(ll0)>sqrt(.Machine$double.eps)) %in%
                     which(abs(l1)>0))
      if(!c0.ok) {
        stop("Provided base correlation is not in the graphpcor model class!")
      }
      x <- ll0[lower.tri(ll0) & (!z0)]
    } else {
      stopifnot(length(x) == nEdges)
      C0 <- vcov(func, theta = x)
    }
    itheta <- which(lower.tri(Q0) & (!is.zero(Q0)))
    ## hessian uses graphpcor:::KLD10
    H <- Hcorrel(
      theta = x,
      p = n,
      parametrization = "itp",
      itheta = itheta,
      d0 = n:1,
      C0 = C0,
      decomposition = decomposition)
    return(H)
  }
#' @describeIn graphpcor
#' The `cgeneric` method for `graphpcor` uses [cgeneric_graphpcor()]
#' @importFrom INLAtools cgeneric
#' @export
cgeneric.graphpcor <- function(model, ...) {
  args <- list(...)
  args$model <- model
  do.call(what = 'cgeneric_graphpcor',
          args = args)
}
#' @describeIn graphpcor
#' The `cgeneric` method for `matrix` uses [cgeneric_graphpcor()]
#' @export
cgeneric.matrix <- function(model, ...) {
  args <- list(...)
  args$model <- graphpcor(model)
  do.call(what = 'cgeneric_graphpcor',
          args = args)
}
