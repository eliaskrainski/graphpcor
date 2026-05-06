#' graphpcor: correlation from nodes and edges
#' @description A `graphpcor` is a graph where
#' a node represents a variable and an edge
#' represent a conditional distribution.
#' The correlation built from a `graphpcor` consider
#' the parameters for the Cholesky of a precision matrix,
#' whose non-zero pattern is given from the graph.
#' @param ... matrix or Matrix (treated as binary) or a
#' vector or list of formula (or character interpreted as formula)
#' @export
graphpcor <- function(...) {
  UseMethod("graphpcor")
}
#' @describeIn graphpcor
#' Each element may be a character or formula.
#' @export
graphpcor.list <- function(...) {
  fch <- lapply(unlist(list(...)), function(x)
    gsub("\"", "", deparse(x)))
  do.call(
    what = "graphpcor",
    args = fch
  )
}
#' @describeIn graphpcor
#' Each term represents a node, and each `~` an edge.
#' @export
graphpcor.formula <- function(...) {
  fch <- lapply(list(...), function(x)
    gsub("\"", "", deparse(x)))
  do.call(
    what = "graphpcor",
    args = fch
  )
}
#' @describeIn graphpcor
#' Each term represents a node, and each `~` an edge.
#' @export
graphpcor.character <- function(...) {
  fch <- sapply(unlist(list(...)), function(x)
    gsub("\"", "", deparse(x)))
  m <- length(fch)
  if(m<1)
    stop("Please provide an argument!")
  ch <- lapply(fch, function(x)
    as.character(as.formula(x)))
  ## right side check, and collect terms
  terms.r <- lapply(ch, function(x) {
    k <- length(x)
    x <- gsub(" ", "", x[k])
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
    jj <- setdiff(jj, ii)
    if((length(ii)>0) & (length(jj)>0)) {
      grel[ii, jj] <- 1
      grel[jj, ii] <- 1
    }
  }
  class(fch) <- 'graphpcor'
  attr(fch, 'nodes') <- allNodes
  attr(fch, 'graph') <- Sparse(grel)
  return(fch)
}
#' @describeIn graphpcor
#' Build a `graphpcor` from a matrix object
#' @export
graphpcor.matrix <- function(...) {
  return(graphpcor(Sparse(list(...)[[1]])))
  ## bellow old code
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
  for(i in 1:nj) {
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
#' @export
graphpcor.Matrix <- function(...) {
  a <- list(...)[[1]]
  diag(a) <- 0
  a <- Sparse(a, na.rm = TRUE, zeros.rm = TRUE)
  a@x <- rep(1.0, length(a@x))
  stopifnot(all.equal(a, t(a), check.attributes = FALSE))
  x <- upperPadding(a)
  nNames <- colnames(a)
  if(is.null(nNames))
    nNames <- rownames(a)
  n <- ncol(x)
  if(is.null(nNames)) {
    nNames <- paste0("x", 1:n)
  }
  dimnames(a) <- list(nNames, nNames)
  ch <- vector("character", n)
  ii <- which(x@j>x@i)
  if(length(ii)>0) {
    sj <- split(nNames[x@j[ii]+1],
                factor(nNames[x@i[ii]+1], nNames))
    for(i in 1:n) {
      nj <- length(sj[[i]])
      if(nj==0) {
        ch[i] <- paste("~", nNames[i])
      } else {
        if(nj==1) {
          ch[i] <- paste(nNames[i], "~", nNames[sj[[i]]])
        } else {
          ch[i] <- paste(nNames[i], "~",
                         paste(nNames[sj[[i]]], collapse = "+"))
        }
      }
    }
  } else {
    ch <- paste("~", nNames)
  }
  class(ch) <- "graphpcor"
  attr(ch, "nodes") <- nNames
  attr(ch, "graph") <- a
  return(ch)
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
#' @export
plot.graphpcor <- function(x, y, ...) {
    ne <- dim(x)
    nodes <- attr(x, "nodes")
    stopifnot(!is.null(nodes))
    stopifnot(ne[1]==length(nodes))
    edgl <- edges(x)
    dotArgs <- list(...)
    if(is.null(dotArgs$Rgraphviz) ||
       !dotArgs$Rgraphviz) { ## depends on igraph
      higraph <- try(do.call(
        what = "require",
        args = list(package = "igraph")), silent = TRUE)
      if(inherits(higraph, "try-error")) {
        cat(higraph)
        stop("Please install 'igraph'!")
      }
      a <- upperPadding(attr(x, "graph"))
      g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = a
      )
      if(is.null(list(...)$arrow.mode)) {
        plot(g, edge.arrow.mode = 0, ...)
      } else {
        plot(g, ...)
      }
    } else {
      havegraph <- do.call(
        what = "require",
        args = list(package = "Rgraphviz"))
      if(havegraph) {
        gr <- do.call(
          what = "graphNEL",
          args = list(nodes = nodes,
                      edgeL = edgl,
                      edgemode = 'undirected')
        )
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

        ag <- do.call(
          what = "agopen",
          args = list(
            graph = gr,
            name = "",
            nodeAttrs = nattr)
        )
        for(k in 1:length(ag@AgEdge)) {
          ag@AgEdge[[k]]@color <- "red"
        }
        getMethod("plot", "Ragraph")(ag)
      } else {
        stop("'Rgraphviz' is not available!")
      }
    }
  }
#' @describeIn graphpcor
#' The `vcov` method for a `graphpcor`
#' @export
vcov.graphpcor <-
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
    iLtheta <- which(lower.tri(G) & (!is.zero(G)))
    LQ0 <- Lprec0(
      theta = theta[-(1:ne[1])],
      p = ne[1],
      iLtheta = iLtheta,
      d0 = d0)

    ## std
    V <- chol2inv(t(LQ0))
    si <- exp(theta[1:ne[1]]) / sqrt(diag(V))
    V <- diag(si) %*% V %*% diag(si)
    dimnames(V) <- names2

    return(V)
}
#' @describeIn graphpcor
#' The precision method for 'graphpcor'
#' @param model graphpcor model object
#' @export
gphcQ <- function(model, ...) {
  V <- vcov(model, ...)
  Q <- chol2inv(chol(V))
  Q[is.zero(Q)] <- 0
  return(Sparse(forceSymmetric(Matrix(Q)),
                zeros.rm = TRUE))
}
#' @describeIn graphpcor
#' The `cgeneric` method for `graphpcor` uses [cgeneric_graphpcor()]
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
