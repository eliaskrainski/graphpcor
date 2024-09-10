#' Set as a graph with nodes representing two kind of
#' variables: children and parent. See details.
#' @param ... a list of formula used as relationship
#' to define the Directed Tree Graph - DTG.
#' Parent nodes shall be in the right side while children
#' (or parents with a parent) in the left side.
#' Please see the examples.
#' @details
#' The children variables are those with an ancestor (parent),
#' and are identified as `c1`, ..., `cn`, where `n` is the
#' total number of children variables.
#' The variables are identified as `p1`, ..., `pm`,
#' where the `m` is the number of parent variables.
#' The main parent (fist) should be identified as `p1`.
#' Parent variables (exept `p1`) have an ancestor,
#' which is a parent variable.
#' @export
#' @examples
#' g1 <- dtg(p1 ~ c1 + c2 - c3)
#' g1
#' summary(g1)
#' plot(g1)
#' precision(g1)
#' precision(g1, theta = 0)
#'
#' g2 <- dtg(p1 ~ c1 + c2 + p2,
#'           p2 ~ c3 - c4)
#' g2
#' summary(g2)
#' plot(g2)
#' precision(g2)
#' precision(g2, theta = c(0, 0))
#'
#' g3 <- dtg(p1 ~ -p2 + c1 + c2,
#'           p2 ~ c3)
#' g3
#' summary(g3)
#' plot(g3)
#' precision(g3)
#' precision(g3, theta = c(0,0))
dtg <- function(...) {

  fch <- as.character(match.call())[-1]
##  fch <- lapply(fch0, as.formula)
  if(length(fch)<1)
    stop("Please provide an argument!")

  ch <- lapply(fch, function(x)
    as.character(as.formula(x)))

  ## left side check
  pp <- sapply(ch, function(x) x[2])
  if (length(pp)>length(unique(pp))) {
    stop("Non-unique parent definition!")
  }
  if(!all(sapply(pp, substr, 1, 1) == "p"))
    stop("Please use the letter 'p' for parents!")
  if(any(is.na(as.integer(
    sapply(ch, function(x)
      substring(x[2], 2)))))) {
    stop("Please use integer after letter 'p' for parents!")
  }

  ## left side check

  ## order
  opp <- order(pp)
  fch <- fch[opp]
  ch <- ch[opp]
  pp <- pp[opp]
  m <- length(pp)
  P <- sort(pp)

  ## right side check, and collect terms
  terms.l <- vector("list", m)
  terms.i <- vector("list", m)
  terms.s <- vector("list", m)
  for(i in 1:m) {
    x <- gsub(" ", "", ch[[i]][3])
##    print(c(x = x))
    if(substr(x, 1, 1) != "-")
      x <- paste0("+", x)
  ##  print(c(x = x))
    i0 <- gregexpr("-", x, fixed = TRUE)[[1]]
    i1 <- gregexpr("+", x, fixed = TRUE)[[1]]
    si <- integer(max(i0, i1))
    if(any(i0>0))
      si[i0] <- -1
    if(any(i1>0))
      si[i1] <- +1
    ##print(si)
    terms.s[[i]] <- si[si != 0]
    schi <- strsplit(x, "-", fixed = TRUE)[[1]]
##    print(schi)
    schi <- unlist(strsplit(schi, "+", fixed = TRUE))
  ##  print(schi)
    if(schi[1]=="") schi <- schi[-1]
    ch.l <- substr(schi, 1, 1)
    ##print(c(ch.l))
    if(!all(ch.l %in% c("p", "c")))
      stop("Invalid variable labeling in ~ ", ch[[i]][3])
    jj.i <- as.integer(substring(schi, 2))
    ##print(jj.i)
    if(any(is.na(jj.i)) | any(jj.i<1))
      stop("Invalid variable numbering in ~ ", ch[[i]][3])
    if(any(ch.l == "p")) {
      if(any(jj.i[ch.l == "p"] <= 1))
        stop("Parent id in ~ ", ch[[i]][3], " must be >", 1)
      if(any(jj.i[ch.l == "p"] > m))
        stop("Wrong parent id in ~ ", ch[[i]][3])
    }
    terms.l[[i]] <- ch.l
    terms.i[[i]] <- jj.i
  }

  rterms <- paste0(
    unlist(terms.l),
    unlist(terms.i))
  urterms <- unique(rterms)
  if(length(rterms) > length(urterms)) {
    stop("Duplicated defintion for ",
         unique(rterms[duplicated(rterms)]))
  }

  ic <- unlist(terms.l) == "c"
  if(any(ic)) {
    n <- max(unlist(terms.i)[ic])
  } else{
    n <- 0L
  }

  trm <- matrix(0L, n+m, m)
  colnames(trm) <- paste0("p", 1:m)
  if(n>0) {
    rownames(trm) <- c(
      paste0("c", 1:n),
      paste0("p", 1:m)
    )
  } else {
    stop("There is no children variable!")
  }

  for(i in 1:m) {
    icc <- terms.l[[i]] == "c"
    if(any(icc)) {
      trm[terms.i[[i]][icc], i] <- terms.s[[i]][icc]
    }
    npc <- sum(!icc)
    if(npc>0) {
      trm[n + terms.i[[i]][!icc], i] <- terms.s[[i]][!icc]
    }
  }

  trm <- trm[-(n+1), , drop = FALSE]

  nr <- rowSums(trm != 0)
  if(any(nr==0)) {
    stop("Missing definition for ",
         rownames(trm)[nr==0])
  }

  class(fch) <- "dtg"
  attr(fch, "childrens") <- n
  attr(fch, "parents") <- m
  attr(fch, "relations") <- trm

  return(fch)

}

#' @export
print.dtg <- function(x, ...) {
  cat("DTG for",
      attr(x, "childrens"),
      "childrens with",
      attr(x, "parents"),
      "parents\n")
  for(i in 1:length(x)) {
    cat(x[[i]], "\n")
  }
}

#' @export
summary.dtg <- function(object, ...) {
  attr(object, "relations")
}

#' @export
dim.dtg <- function(x, ...) {
  trm <- attr(x, "relations")
  m <- ncol(trm)
  c(nrow(trm) - m + 1, m)
}

#' @export
setMethod(
  "drop",
  "dtg",
  function(x) {
    stopifnot((m <- length(x))>1)
    trm0 <- attr(x, "relations")
    stopifnot(ncol(trm0) == m)
    ilast <- which(rownames(trm0) == (colnames(trm0)[m]))
    trm <- trm0[-ilast, 1:(m-1), drop = FALSE]

    args <- lapply(1:(m-1), function(i) {
      j <- trm[, i] !=0
      s <- ifelse(trm[j, i] < 0, "-", "+")
      if(s[1] == "+") s[1] <- ""
      e <- paste(
        colnames(trm)[i],
        "~",
        paste(s, rownames(trm)[j], collapse = "")
      )
      update(y ~ x, e)
    }
    )
    do.call(
      "dtg",
      args)
  }
)

#' @export
setMethod(
  "edges",
  "dtg",
  function(object, which, ...) {

    trm <- attr(object, "relations")
    m <- ncol(trm)
    n <- nrow(trm)-m+1

    edgl <- vector("list", m + n)
    names(edgl) <- c(paste0("p", 1:m),
                     paste0("c", 1:n))
    for(i in 1:m) {
      w1 <- trm[, i] != 0
      edgl[[i]] <- list(
        n = sum(w1),
        edges = rownames(trm)[w1],
        weights = trm[w1, i]
      )
      edgl[[i]]$term <- edgl$edges
      edgl[[i]]$parent <- substr(edgl[[i]]$edges, 1, 1) == "p"
      edgl[[i]]$id <- as.integer(substring(edgl[[i]]$edges, 2))
      edgl[[i]]$signal <- edgl[[i]]$weights
    }
    return(edgl)
  }
)

#' @export
setMethod(
  "plot",
  "dtg",
  function(x, y, ...) {

    ## the graph
    edgl <- edges(x)
    nodes <- names(edgl)
    gr <- graph::graphNEL(
      nodes = nodes,
      edgeL = edgl,
      edgemode='directed')

    trm <- attr(x, "relations")
    m <- ncol(trm)
    n <- nrow(trm)-m+1

    mc <- lapply(
      match.call(expand.dots = TRUE)[-1],
      eval)
    nargs <- names(mc)

    nattr <- lapply(
      list(
        color =  {
          if(any(nargs == "color"))
            mc$color[1:2]
          else
            c("red", "blue")
          },
        fillcolor = {
          if(any(nargs == "fillcolor"))
            mc$fillcolor[1:2]
          else
            c("lightsalmon", "lightblue")
          },
        shape = {
          if(any(nargs == "shape"))
             mc$shape[1:2]
           else
             c("box", "circle")
        },
        height = {
          if(any(nargs == "height"))
            mc$height[1:2]
          else
            c(0.5, 0.5)
        },
        width = {
          if(any(nargs == "width"))
            mc$width[1:2]
          else
            c(0.75, 0.75)
        },
        fontsize = {
          if(any(nargs == "fontsize"))
            mc$fontsize
          else
            c(14, 14)
        }
      ), rep, times = c(m, n))
    for(i in 1:length(nattr))
      names(nattr[[i]]) <- nodes

    # eattr <- buildEdgeList(gr)
    # for(i in 1:m) {
    #   j <- which(trm[, i] != 0)
    #   for(k in 1:length(j)) {
    #     nam <- paste0(colnames(trm), "~", rownames(trm)[j[k]])
    #     l <- which(names(eattr) == nam)
    #     eattr[[l]]@attrs$weight <- paste(trm[j[k], i])
    #   }
    # }
    plot(gr,
         nodeAttrs = nattr,
         ##edgeAttrs = eattr,
         ...)

  }
)

#' @export
variance <- function(x, ...) {
  UseMethod("variance")
}

#' @export
variance.default <- function(x, ...) {
  return(var(x))
}

#' @export
precision <- function(x, ...) {
  UseMethod("precision")
}

#' @export
precision.default <- function(x, ...) {
  v <- variance(x, ...)
  return(
    forwardsolve(
      backsolve(
        chol(v)
        )
      )
    )
}

#' @export
precision.dtg <- function(x, ...) {
  edgl <- edges(x)
  q.el <- edtg2precision(
    edgl[1:length(x)])
  mc <- list(...)
  nargs <- names(mc)
  Q <- q.el$q
  if(any(nargs == "theta")) {
    nc <- q.el$nc
    Q[q.el$iq2th] <- Q[q.el$iq2th] +
      exp(-2 * mc$theta[q.el$i2th])
    Q[q.el$iq1th] <- -1.0 * q.el$sth *
      exp(-mc$theta[q.el$i1th])
  }
  d <- dim(x)
  rownames(Q) <- colnames(Q) <-
    names(edgl)[c(d[2] + 1:d[1], 1:d[2])]
  return(Q)
}

#' Internal function to extract elements to
#' build the precision from the DTG edges.
#' @param d.el the list of the first n edges of a DTG.
edtg2precision <- function(d.el) {
  stopifnot(all(substr(names(d.el), 1, 1) == "p"))
  stopifnot(length(d.el) == length(unique(names(d.el))))
  ip <- as.integer(substring(names(d.el), 2))
  np <- length(ip)
  stopifnot(np == length(unique(ip)))
  nc <- sum(sapply(d.el, function(x) sum(!x$parent)))
  p.nc <- sapply(d.el, function(x) x$n)
  dd <- c(rep(1, nc), p.nc)
  stopifnot((nc + np) == length(dd))
  q0 <- diag(x = dd, nrow = nc + np, ncol = nc + np)
  ij <- matrix(1:((nc+np)^2), nc+np, nc+np)
  iq1th <- integer(2 * (np - 1))
  sch <- iq1ch <- integer(2*nc)
  sth <- i1th <- integer(np-1)
  iq2th <- i2th <- integer(np)
  k0 <- 0
  k2 <- k1 <- 0
  for(i in 1:np) {
    i0 <- which(!d.el[[i]]$parent)
    nci <- length(i0)
    if(nci>0) {
      j <- d.el[[i]]$id[i0]
      sch[k0 + 1:nci] <- d.el[[i]]$signal[i0]
      iq1ch[k0 + 1:nci] <- ij[(col(ij) == (nc+i)) & (row(ij) %in% j)]
      k0 <- k0 + nci
      sch[k0 + 1:nci] <- d.el[[i]]$signal[i0]
      iq1ch[k0 + 1:nci] <- ij[(row(ij) == (nc+i)) & (col(ij) %in% j)]
      k0 <- k0 + nci
      q0[j, nc+i] <- -d.el[[i]]$signal[i0]
      q0[nc+i, j] <- -d.el[[i]]$signal[i0]
    }
    i2th[k1 + 1] <- i
    iq2th[k1 + 1] <- ij[(col(ij) == (nc+i)) & (row(ij) == (nc+i))]
    k1 <- k1 + 1
    i0 <- which(d.el[[i]]$parent)
    nj <- length(i0)
    if(nj>0) {
      j0 <- d.el[[i]]$id[i0]
      i1th[k2 + 1:nj] <- j0
      sth[k2 + 1:nj] <- d.el[[i]]$signal[i0] ## carry on the signal
      j <- nc + j0
      iq1th[k2 + 1:nj] <- ij[, nc+i][j]
      k2 <- k2 + nj
      i1th[k2 + 1:nj] <- j0
      sth[k2 + 1:nj] <- d.el[[i]]$signal[i0] ## carry on the signal
      iq1th[k2 + 1:nj] <- ij[nc+i, ][j]
      k2 <- k2 + nj
      q0[j, nc+i] <- -d.el[[i]]$signal[i0]
      q0[nc+i, j] <- -d.el[[i]]$signal[i0]
    }
  }
  stopifnot(k1 == np)
  stopifnot(k2 == (2*(np-1)))
  return(list(
    nc = as.integer(nc),
    np = as.integer(np),
    i2th = as.integer(i2th),
    iq2th = as.integer(iq2th),
    i1th = as.integer(i1th),
    iq1th = as.integer(iq1th),
    iq1ch = as.integer(iq1ch),
    sch = as.double(sch),
    sth = as.double(sth),
    q = q0
  ))
}

#' @export
variance.dtg <- function(x, ...) {
  mc <- lapply(
    match.call(
      expand.dots = TRUE)[-1],
    eval)
  nargs <- names(mc)
  nm <- dim(x)
  edgl <- edges(x)
  ij <- edtg2variance(edgl[1:nm[2]])
  np <- length(ij$iv)
  nc <- length(ij$iparent)
  stopifnot(all(c(nc, np) == nm))
  if(any(nargs == "theta")) {
    vi <- sapply(ij$iv, function(i)
      sum(exp(2 * mc$theta[i])))
  } else {
    vi <- sapply(ij$iv, function(i) length(i))
  }
  vv <- diag(nc) + vi[ij$itop]
  rownames(vv) <- colnames(vv) <- names(edgl)[np + 1:nc]
  return(t(vv * ij$schildren) * ij$schildren)
}

#' Internal function to extract elements to
#' build the covariance matrix from the DTG edges.
#' @param d.el the list of the first n edges of a DTG.
edtg2variance <- function(d.el) {
  np <- length(d.el)
  iv <- lapply(1:np, function(i) i)
  for(i in 1:np) {
    ip <- which(d.el[[i]]$parent)
    if(length(ip)>0) {
      jj <- d.el[[i]]$id[ip]
      for(j in jj) {
        iv[[j]] <- c(iv[[j]], iv[[i]])
      }
    }
  }
  NC <- sum(sapply(d.el, function(s) sum(!s$parent)))
  sch <- integer(NC)
  iP <- integer(NC)
  for(i in 1:np) {
    ic <- which(!d.el[[i]]$parent)
    if(length(ic)>0) {
      sch[d.el[[i]]$id[ic]] <-
        d.el[[i]]$signal[ic]
      for(j in 1:length(ic)) {
        iP[d.el[[i]]$id[ic[j]]] <- i
      }
    }
  }
  itop <- matrix(0L, NC, NC)
  for(i in 1:NC) {
    for(j in 1:NC) {
      itop[i, j] <- max(intersect(iv[[iP[i]]], iv[[iP[j]]]))
    }
  }
  stopifnot(all.equal(iP,diag(itop)))
  return(list(iparent = iP, iv = iv, itop = itop, schildren=sch))
}

#' cgeneric method in `INLA` from a list of expressions
#' defining a Direct Tree Graph - DTG correlation model
#' to be used as a model in a `INLA` `f()` model component.
#'
#' @param x the DTG model specification.
#' @param sigma.prior.reference a vector with the reference values
#' to define the prior for the standard deviation parameters.
#' @param sigma.prior.probability a vector with the probability values
#' to define the prior for the standard deviation parameters.
#' @param lambda the lambda for the graph correlation prior.
#' @param iprior integer to define with prior is going to be used
#' for the correlation. 1 is for normal for v_i,
#' 2 is for pc-precision and 3 is for the derived PC (see the paper).
#' @param useINLAprecomp logical indicating if
#' it is to be used the shared object within INLA.
#' @details
#'  The correlation prior as in the paper depends on the lambda value.
#'  The prior for each \eqn{sigma_i} is the Penalized-complexity prior
#' which can be defined from the following probability statement
#'  P(sigma > U) = a.
#' where "U" is a reference value and "a" is a probability.
#' The values "U" and probabilities "a" for each \eqn{sigma_i}
#' are passed in the `sigma.prior.reference` and `sigma.prior.probability`
#' arguments.
#' If a=0 then U is taken to be the fixed value of the corresponding sigma.
#' E.g. if there are three sigmas in the model and one supply
#'  sigma.prior.reference = c(1, 2, 3) and
#'  sigma.prior.probability = c(0.05, 0.0, 0.01)
#' then the sigma is fixed to 2 and not estimated.
#' @return objects to be used in the f() formula term in INLA.
#' @export
cgeneric.dtg <-
  function(dtg,
           sigma.prior.reference,
           sigma.prior.probability,
           lambda,
           iprior = 3,
           debug = FALSE,
           useINLAprecomp = !TRUE) {


    dd <- dim(dtg)
    d.el <- edges(dtg)[1:dd[2]]
    ich <- unlist(lapply(d.el, function(x)
      x$id[!x$parent]))
    sch <- unlist(lapply(d.el, function(x)
      x$signal[!x$parent]))
    sch <- sch[ich]
    if(debug) {
      cat(c(sch = sch), "\n")
    }
    d.elc <- edtg2variance(d.el)
    if(debug) {
      print(str(d.elc))
    }
    np <- length(dtg)
    nv <- sapply(d.elc$iv, length)
    if(debug)
      cat("np = ", np, " and nv: ", nv, "\n")
    iiv <- rep(1:np, nv)
    jjv <- unlist(lapply(d.elc$iv, sort))
    itop <- d.elc$itop
    if(debug) {
      cat(c(iiv=iiv), "\n")
      cat(c(jjv=jjv), "\nitop:\n")
      print(itop)
    }
    nc <- nrow(itop)
    ii <- col(itop)[!upper.tri(itop)]
    jj <- row(itop)[!upper.tri(itop)]
    if(debug) {
      print(str(list(nc=nc,ii=ii,jj=jj)))
    }

    stopifnot(length(sigma.prior.reference) == nc)
    stopifnot(length(sigma.prior.probability) == nc)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    stopifnot(lambda>0)
    stopifnot(iprior %in% (1L:3L))

      the_model <- cgeneric.default(
          model = "inla_cgeneric_corgraphs_sfixed",
          debug = as.logical(debug),
          useINLAprecomp = as.logical(useINLAprecomp),
          n = as.integer(nc),
          np = as.integer(np),
          nv = as.integer(nv),
          ipar = as.integer(d.elc$iparent-1L),
          iiv = as.integer(iiv-1L),
          jjv = as.integer(jjv-1L),
          itop = as.integer(itop-1L),
          ii = as.integer(ii-1L),
          jj = as.integer(jj-1L),
          iprior = as.integer(iprior),
          lambda = as.double(lambda),
          slambdas = as.double(slambdas),
          schildren = as.double(sch)
      )

    return(the_model)

  }
