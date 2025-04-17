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
  nodes <- unique(unlist(lapply(1:m, function(i)
    c(nodesL[i], terms.r[[i]]))))
  nNodes <- length(nodes)
  grel <- matrix(
    0, m, length(nodesR),
    dimnames = list(nodesL, nodesR))
  for(i in 1:m) {
    jj <- pmatch(terms.r[[i]], nodesR)
    grel[i, jj] <- 1
  }
  class(fch) <- 'graphpcor'
  attr(fch, 'nodes') <- nodes
  attr(fch, 'relationship') <- grel
  return(fch)
}
#' @describeIn graphpcor
#' Build a `graphpcor` from a matrix
#' @importFrom stats as.formula
#' @export
graphpcor.matrix <- function(...) {
  x <- list(...)[[1]]
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
  return(do.call(what = 'graphpcor',
                 args = lapply(argl, as.formula)))
}
#' @describeIn graphpcor
#' The print method for `graphpcor`
#' @param x graphpcor
#' @export
print.graphpcor <- function(x, ...) {
  cat("A graphpcor for",
      length(attr(x, 'nodes')), "variables",
      "with", sum(attr(x, 'relationship')), "edges.\n")
}
#' @describeIn graphpcor
#' The summary method for `graphpcor`
#' @param object graphpcor
#' @export
summary.graphpcor <- function(object, ...) {
  attr(object, "relationship")
}
#' @describeIn graphpcor
#' The dim method for `graphpcor`
#' @param x graphpcor
#' @export
dim.graphpcor <- function(x, ...) {
  c(nodes=length(attr(x, 'nodes')),
    edges=sum(attr(x, 'relationship')))
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
  grel <- attr(graph, "relationship")
  L <- matrix(
    0.0, ne[1], ne[1],
    dimnames = list(nodes, nodes)
  )
  for(i in 1:nrow(grel)) {
    ii <- pmatch(rownames(grel)[i], nodes)
    jj <- pmatch(colnames(grel)[grel[i, ]!=0], nodes)
    L[ii, jj] <- (-1)
  }
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  return(L)
}
#' @describeIn graphpcor
#' Build the unite diagonal lower triangle matrix
#' @param x a `graphpcor` object
setMethod(
  "chol",
  "graphpcor",
  function(x, ...) {
    ne <- dim(x)
    G <- Laplacian(x)
    idx <- which(lower.tri(G) & (!is.zero(G)))
    L <- diag(x = 1:ne[1], nrow = ne[1], ncol = ne[1])
    stopifnot(length(idx)==ne[2])
    args <- list(...)
    stopifnot(length(args$theta)==ne[2])
    L[idx] <- args$theta
    ll <- t(chol(G + diag(ne[1])))
    ifill <- which(is.zero(G) & (!is.zero(ll)))
    if(length(ifill)>0) {
      L <- fillLprec(L, ifill)
    }
    return(t(L))
  }
)
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
    Q <- Laplacian(object)
    stopifnot(ne[1]==nrow(Q))
    stopifnot((2*ne[2])==(sum(!is.zero(Q))-ne[1]))
    if(length(theta)==ne[2]) {
      theta <- c(rep(0.0, ne[1]), theta)
    } else {
      stopifnot(length(theta)==sum(ne))
    }
    L <- getMethod('chol', 'graphpcor')(
      object, theta = theta[-(1:ne[1])])
    V <- chol2inv(L)
    si <- exp(theta[1:ne[1]]) / sqrt(diag(V))
    V <- diag(si) %*% V %*% diag(si)
    return(V)
  }
)
#' @describeIn graphpcor
#' The precision method for 'graphpcor'
#' @param model graphpcor model object
#' @export
prec.graphpcor <- function(model, ...) {
  ne <- dim(model)
  Q <- Laplacian(model)
  stopifnot(ne[1]==nrow(Q))
  stopifnot((2*ne[2])==(sum(!is.zero(Q))-ne[1]))
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    return(chol2inv(chol(
      vcov(model, ...))
    ))
  } else {
    warning("missing `theta`, returning Laplacian!")
  }
  return(Q)
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
#' @export
hessian.graphpcor <-
  function(func,
           x,
           method = "Richardson",
           method.args = list(),
           ...) {
    decomposition <- c("eigen", "svd", "chol")
    if(is.null(list(...)$decomposition)) {
      decomposition <- "eigen"
    } else {
      decomposition <- match.arg(
        list(...)$decomposition,
        decomposition)
    }
    Q0 <- Laplacian(func)
    nEdges <- sum((!is.zero(Q0)) & lower.tri(Q0))
    z0 <- is.zero(Q0)
    n <- nrow(Q0)
    l1 <- t(chol(Q0 + diag(1.0, n, n)))
    if(inherits(x, "matrix")) {
      ## maybe optim() to get theta that
      ## give it close to C0?
      ## For now check the elements of L from C0^{-1}
      C0 <- cov2cor(x)
      qq0 <- chol2inv(chol(C0))
      ll0 <- t(chol(qq0))
      c0.ok <- all(which(abs(ll0)>sqrt(.Machine$double.eps)) %in%
                     which(abs(l1)>0))
      if(!c0.ok) {
        stop("Provided base correlation not in the graphpcor model class!")
      }
      for(i in 1:nrow(C0))
        ll0[i, ] <- ll0[i, ]/ll0[i, i]
      x <- ll0[lower.tri(ll0) & (!z0)]
    } else {
      stopifnot(length(x) == nEdges)
      C0 <- vcov(func, theta = x)
    }
    ## hessian uses graphpcor:::KLD10
    H <- hessian(
      func = function(th)
        KLD10(vcov(func, theta = th), C0),
      x = x,
      method = method,
      method.args = method.args,
      ...)
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
      hn.5 <- chol(hn, pivot = TRUE)
      hneg.5 <- matrix(hn.5[, order(attr(hn.5, "pivot"))], nrow(H))
    }
    stopifnot(all.equal(H, tcrossprod(h.5)))
    attr(H, "base") <- x
    attr(H, "h.5") < h.5
    attr(H, "hneg.5") <- hneg.5
    attr(H, "decomposition") <- Hd
    return(H)
  }
#' @describeIn cgeneric
#' The `cgeneric` method for `graphpcor` uses [cgeneric_graphpcor()]
#' @export
cgeneric.graphpcor <- function(...) {
  args <- list(...)
  args$graph <- args$model
  args$model <- NULL
  do.call(what = 'cgeneric_graphpcor',
          args = args)
}
