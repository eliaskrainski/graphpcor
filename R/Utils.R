#' Function to remove the last parent from the given DAG
#' @param dag model structure given as a formula list
#' @return a dag without the last parent whose children
#' as children of its parent.
#' @export
#' @examples
#' dag1 <- list(
#'      p1 ~ p2 - p3 + c1 + c2,
#'      p2 ~ c3 - c4,
#'      p3 ~ -c5 + c6)
#' dag2 <- dag_rlp(dag1)
#' dag2
#' dag_rlp(dag2)
dag_rlp <- function(dag) {
  np <- length(dag)
  stopifnot(np>1)
  stilde <- sapply(dag, function(x)
    strsplit(gsub(" ", "", as.character(x)),
             split = "~"))
  p.rm <- stilde[2, np]
  il.rm <- grep(p.rm, stilde[3, ])
  ii <- gregexpr(p.rm, stilde[3, il.rm])[[1]]
  iin <- ii+attr(ii, "match.length")
  if(ii>2) {
    l.up0 <- substr(stilde[3, il.rm], 1, ii-2)
    if(iin<nchar(stilde[3, il.rm]))
      l.up0 <- paste0(l.up0,
                      substring(stilde[3, il.rm], iin))
  } else {
    l.up0 <- substring(stilde[3, il.rm], iin+1)
  }
  if(substr(stilde[3, np], 1, 1) != "-") {
     stilde[3, il.rm] <- paste(l.up0, "+", stilde[3, np])
  } else {
    stilde[3, il.rm] <- paste(l.up0, stilde[3, np])
  }
  il.env <- attr(dag[[il.rm]], ".Environment")
  dag[[il.rm]] <- as.formula(
    paste(stilde[2, il.rm], "~", stilde[3, il.rm]))
  attr(dag[[il.rm]], ".Environment") <- il.env
  return(dag[1:(np-1)])
}
#' Function to extract the elements of a DAG tree
#' @return a list with same length of S.
#' Each element as a named list to specify the left hand side
#' of the graph terms.
#' \itemize{
#'  \item parent: a logical vector indicating if the term is parent
#'  \item id: the number identifying the term
#'  \item signal: either minus or plus one
#'  }
dag_elements <- function(dag, debug = FALSE) {
  np <- length(dag)
  stilde <- sapply(dag, function(x)
    strsplit(gsub(" ", "", as.character(x)),
             split = "~"))

  if(debug>1)
    print(stilde)
  idParents <- unlist(stilde[2, ])
  stopifnot(length(dag) == length(unique(idParents)))
  stopifnot(all(substr(idParents, 1, 1) == "p"))

  iParents1 <- as.integer(substring(idParents, 2))
  if(debug)
    print(iParents1)
  NP <- length(iParents1)
  if(debug)
    cat("NP = ", NP, "\n")

  stilde

  dag_sl <- vector("list", np)
  names(dag_sl) <- idParents
  for(k in 1:np) {
    if(debug>1)
      cat("k =", k, "\n")
    ch <- stilde[3, k]
    if(substr(ch, 1, 1) != "-")
      ch <- paste0("+", ch)
    p <- s <- integer(nchar(ch))
    i1 <- gregexpr("-", ch, fixed = TRUE)[[1]]
    if(debug>1)
      print(new("integer", i1))
    if(any(i1>0)) {
      s[i1[i1>0]] <- -1
      p[i1[i1>0]] <- i1[i1>0]
    }
    i2 <- gregexpr("+", ch, fixed = TRUE)[[1]]
    if(debug>1)
      print(new("integer", i2))
    if(any(i2>0)) {
      s[i2[i2>0]] <- 1
      p[i2[i2>0]] <- i2[i2>0]
    }
    ss <- strsplit(substring(ch, 2),
                   split = "+", fixed = TRUE)[[1]]
    ss <- unlist(lapply(ss, function(x)
      strsplit(x, "-", fixed = TRUE)[[1]]))
    if(debug>1)
      print(ss)
    dag_sl[[k]] <- list(
      n = length(ss),
      term = ss,
      parent = substr(ss, 1, 1) == "p",
      id = as.integer(substring(ss, 2)),
      signal = s[s!=0]
    )
    if(debug)
      print(str(dag_sl[[k]]))
  }
  return(dag_sl)
}
#' Function to build the precision elements from
#' the elements extracted from the DAG tree.
#' @return the elements to build the precision matrix.
dag_e2precision <- function(d.el) {
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
#' Function to build the precision matrix for DAG tree.
#' @param theta log of the conditional variances of the
#' parent variables
#' @param s.children vector s (one for each children).
#' Default is null and it is taken as (+1 or -1) from the formula.
#' @return a precision matrix for given DAG tree and parameters
#' @export
#' @examples
#' dag1 <- list(
#'      p1 ~ p2 - p3 + c1 + c2,
#'      p2 ~ c3 + c4,
#'      p3 ~ -c5 + c6)
#' Q1 <- dag_precision(dag1, theta = c(0, 0, 0))
#' cov1 <- chol2inv(chol(Q1))
#' round(100 * cov2cor(cov1[1:6, 1:6]))
#'
#' dag2 <- list(
#'      p1 ~ p2 + c1 + c2,
#'      p2 ~ -p3 + c3 + c4,
#'      p3 ~ -c5 + c6)
#' Q2 <- dag_precision(dag2, theta = c(0, 0, 0))
#' cov2 <- chol2inv(chol(Q2))
#' round(100 * cov2cor(cov2[1:6, 1:6]))
dag_precision <- function(dag, theta, s.children = NULL, debug = FALSE, new = TRUE) {
  if(new) {
    d.el <- dag_elements(dag)
    q.el <- dag_e2precision(d.el)
    Q <- q.el$q
    nc <- q.el$nc
    Q[q.el$iq2th] <- Q[q.el$iq2th] +
        exp(-2 * theta[q.el$i2th])
    Q[q.el$iq1th] <- -1.0 * q.el$sth * exp(-2*theta[q.el$i1th])
  } else {
    np <- length(dag)
    stilde <- sapply(dag, function(x)
      strsplit(gsub(" ", "", as.character(x)),
               split = "~"))
    if(debug)
      print(stilde)
    stopifnot(all(substr(stilde[2, ], 1, 1) == "p"))

    iParents1 <- as.integer(substring(
      unlist(stilde[2, ]), 2))
    if(debug)
      print(iParents1)
    NP <- length(iParents1)
    if(debug)
      cat("NP = ", NP, "\n")

    dag_sl <- dag_elements(dag, debug = debug)
    NC <- sum(sapply(dag_sl, function(s)
      sum(!s$parent)))
    if(debug)
      cat("NC = ", NC, "\n")

    q <- exp(theta)

    Qa <- cbind(
      rbind(diag(NC), matrix(0, NP, NC)),
      rbind(matrix(0, NC, NP), diag(q)))
    Qth <- cbind(
      rbind(diag(NC), matrix(0, NP, NC)),
      rbind(matrix(0, NC, NP), diag(NP)*0.0))
    jj <- NC + iParents1
    for(k in 1:np) {
      Sk <- dag_sl[[k]]
      ich <- !Sk$parent
      if (any(ich)) {
        ii <- Sk$id[ich]
        Qa[jj[k], jj[k]] <- Qa[jj[k], jj[k]] + length(ii)
        Qa[ii, jj[k]] <- -Sk$signal[ich]
        Qa[jj[k], jj[k]] <- Qa[jj[k], jj[k]] + length(ii)
        Qa[ii, jj[k]] <- -Sk$signal[ich]
      }
      if (any(Sk$parent)) {
        i1 <- NC + iParents1[k]
        ii <- NC + Sk$id[Sk$parent]
        sq <- -Sk$signal[Sk$parent] * q[ii-NC]
        for(i in 1:length(ii)) {
          Qa[i1, i1] <- Qa[i1, i1] + q[ii[i] - NC]
          if(ii[i]>i1) {
            Qa[i1, ii[i]] <- Qa[i1, ii[i]] + sq[i]
          } else {
            Qa[ii[i], i1] <- Qa[ii[i], i1] + sq[i]
          }
        }
      }
    }
    if(debug>1)
      print(Qa)
    Q <- Qa
    Q[lower.tri(Q)] <- t(Qa)[lower.tri(Q)]

  }
  return(Q)
}
#' Internal function to extract elements to build
#' the covariance matrix
dag_e2covariance <- function(d.el) {
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
#   return(iv)
   NC <- sum(sapply(d.el, function(s) sum(!s$parent)))
   iP <- integer(NC)
   for(i in 1:np) {
     ic <- which(!d.el[[i]]$parent)
     if(length(ic)>0)
       for(j in 1:length(ic)) {
         iP[d.el[[i]]$id[ic[j]]] <- i
       }
   }
   itop <- matrix(0L, NC, NC)
   for(i in 1:NC) {
     for(j in 1:NC) {
       itop[i, j] <- max(intersect(iv[[iP[i]]], iv[[iP[j]]]))
     }
   }
   stopifnot(all.equal(iP,diag(itop)))
   return(list(iv = iv, iparent = iP, itop = itop))
}
#' Function to compute the covariance between the
#' children variables for a given dag and the
#' parameters of the parents
#' @export
#' @examples
#' dag1 <- list(
#'      p1 ~ p2 + c1 + c2,
#'      p2 ~ -c3 + c4)
#' dag_covariance(dag1, c(0, 0))
#' dag_covariance(dag_rlp(dag1), 0)
dag_covariance <- function(dag, theta, s.children = NULL) {
  stopifnot(length(theta) == length(dag))
  d.el <- dag_elements(dag)
  ij <- dag_e2covariance(d.el)
  np <- length(ij$iv)
  nc <- length(ij$iparent)
  vv <- diag(nc)
  for(i in 1:nc) {
    for (j in 1:nc) {
      vi <- sum(exp(2 * theta[ij$iv[[ij$iparent[i]]]]))
      vj <- sum(exp(2 * theta[ij$iv[[ij$iparent[j]]]]))
      vv[i, j] <- vv[i, j] + vi * vj
    }
  }
  if(!is.null(s.children)) {
    stopifnot(length(s.children) == nc)
    sc <- diag(s.children)
  } else {
    ich <- unlist(lapply(d.el, function(x)
      x$id[!x$parent]))
    sch <- unlist(lapply(d.el, function(x)
      x$signal[!x$parent]))
    sc <- diag(sch[ich])
  }
  return(sc%*%vv%*%sc)
}
