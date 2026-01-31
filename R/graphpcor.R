#' @describeIn graphpcor
#' Each term to represent a node, and
#' each `~` to represent an edge.
#' @param ... a list of arguments
#' @importFrom stats as.formula
#' @export
#' @example demo/graphpcor.R
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
#' Build a `graphpcor` from a matrix object
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
    vnams <- paste0("x", 1:ne[1])
  }
  argl <- list()
  adde <- 0
  nj <- integer(ne[1]-1)
  for(i in 1:(ne[1]-1)) {
    jj <- intersect(
      (i+1):ne[1],
      which(!iz[i, ]))
    nj[i] <- length(jj)
    if(nj[i]>0) {
      argl[[i]] <- paste(
        vnams[i], "~",
        paste(vnams[jj], collapse = " + "))
      adde <- adde + length(jj)
    }
  }
  stopifnot(adde==ne[2])
  argl <- argl[which(nj>0)]
  return(do.call(what = 'graphpcor',
                 args = lapply(argl, as.formula)))
}
#' @describeIn graphpcor
#' Build a `graphpcor` for a Matrix object
#' @importFrom stats as.formula
#' @importFrom INLAtools Sparse upperPadding
#' @export
graphpcor.Matrix <- function(...) {
  x <- upperPadding(Sparse(list(...)[[1]]))
  ne <- c(nrow(x), length(x@x))
  sij <- split(x@j+1L, x@i+1L)
  vnams <- rownames(x)
  if(is.null(vnams)) {
    vnams <- paste0("x", 1:ne[1])
  }
  argl <- list()
  for(i in 1:length(sij)) {
    xi <- vnams[as.integer(names(sij)[i])]
    xj <- setdiff(vnams[sij[[i]]], xi)
    argl[[i]] <- paste(xi, "~", paste(xj, collapse = "+"))
  }
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
#' @export
dim.graphpcor <- function(x, ...) {
  c(nodes=length(attr(x, 'nodes')),
    edges=sum(attr(x, 'graph'))/2)
}
#' @describeIn graphpcor
#' The `plot` method for a `graphpcor`
#' @param y not used
#' @method plot graphpcor
#' @importFrom methods getMethod
#' @export
plot.graphpcor <- function(x, y, ...) {
    ne <- dim(x)
    nodes <- attr(x, "nodes")
    stopifnot(!is.null(nodes))
    stopifnot(ne[1]==length(nodes))
    edgl <- edges(x)
    if(requireNamespace("igraph")) {
      g <- igraph::graph_from_adjacency_matrix(attr(x, "graph"))
      plot(g, ...)
    } else {
      if(requireNamespace("graph")) {
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
    }
  }
#' @describeIn graphpcor
#' The `vcov` method for a `graphpcor`
#' @importFrom methods getMethod
#' @export
setMethod(
  "vcov",
  "graphpcor",
  function(object, ...) {
    ne <- dim(object)
    p <- ne[1]
    m <- ne[2]
    stopifnot(p>0)
    stopifnot(m>0)

    G <- Laplacian(object)
    names2 <- dimnames(G)
    stopifnot(ne[1]==nrow(G))
    stopifnot((2*ne[2])==(sum(!is.zero(G))-ne[1]))

    ## collect theta
    mc <- list(...)
    nargs <- names(mc)
    if(!any(nargs == "theta")) {
      stop("Please provide 'theta'!")
    }
    theta <- mc$theta
    ## setup full theta
    if(length(theta)==ne[2]) {
      theta <- c(rep(0.0, ne[1]), theta)
    } else {
      stopifnot(length(theta)==sum(ne))
    }

    ## collect/define 'd0'
    if(any(nargs == "d0")) {
      d0 <- mc$d0
    } else {
      d0 <- ne[1]:1
    }

    ## build lower Cholesky of Q0
    itheta <- which(lower.tri(G) & (!is.zero(G)))
    LQ0 <- Lprec0(
      theta = theta[-(1:ne[1])],
      p = ne[1],
      itheta = itheta,
      d0 = ne[1]:1)

    ## std
    V <- chol2inv(t(LQ0))
    si <- exp(theta[1:ne[1]]) / sqrt(diag(V))
    V <- diag(si) %*% V %*% diag(si)
    dimnames(V) <- names2

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
  V <- vcov(model, ...)
  Q <- chol2inv(chol(V))
  Q[is.zero(Q)] <- 0
  return(Sparse(forceSymmetric(Matrix(Q)),
                zeros.rm = TRUE))
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
#' The `cgeneric` method for `Matrix` uses [cgeneric_graphpcor()]
#' @export
cgeneric.Matrix <- function(model, ...) {
  args <- list(...)
  args$model <- graphpcor(model)
  do.call(what = 'cgeneric_graphpcor',
          args = args)
}
