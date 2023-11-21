#' Function that returns density for a given model structure
#' @param dag model structure given as a formula list
#' @return a list with same length of S.
#' Each element as a named list to specify the left hand side
#' of the graph terms.
#' \itemize{
#'  \item parent: a logical vector indicating if the term is parent
#'  \item id: the number identifying the term
#'  \item signal: either minus or plus one
#'  }
dag_elements <- function(dag, debug = FALSE) {
  nS <- length(dag)
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

  dag_sl <- vector("list", nS)
  names(dag_sl) <- idParents
  for(k in 1:nS) {
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
#' Function to build the precision elements
#' @param dag model structure given as a formula list
dag_precision_elements <- function(dag, UPLO = FALSE) {
  el <- dag_elements(dag)
  stopifnot(all(substr(names(el), 1, 1) == "p"))
  ip <- as.integer(substring(names(el), 2))
  stopifnot(length(ip) == length(unique(ip)))
  p <- length(ip)
  n <- sum(sapply(el, function(x) sum(!x$parent)))
  p.nc <- sapply(el, function(x) x$n)
  dd <- c(rep(1, n), p.nc)
  stopifnot((n+p) == length(dd))
  q0 <- diag(x = dd, nrow = n + p, ncol = n + p)
  ij <- matrix(1:((n+p)^2), n+p, n+p)
  iq1th <- integer(2 * (p - 1))
  sch <- iq1ch <- integer(2*n)
  sth <- i1th <- integer(p-1)
  iq2th <- i2th <- integer(p)
  k0 <- 0
  k2 <- k1 <- 0
  for(i in 1:p) {
    i0 <- which(!el[[i]]$parent)
    nci <- length(i0)
    if(nci>0) {
      j <- el[[i]]$id[i0]
      sch[k0 + 1:nci] <- el[[i]]$signal[i0]
      iq1ch[k0 + 1:nci] <- ij[(col(ij) == (n+i)) & (row(ij) %in% j)]
      k0 <- k0 + nci
      sch[k0 + 1:nci] <- el[[i]]$signal[i0]
      iq1ch[k0 + 1:nci] <- ij[(row(ij) == (n+i)) & (col(ij) %in% j)]
      k0 <- k0 + nci
      q0[j, n+i] <- -el[[i]]$signal[i0]
      q0[n+i, j] <- -el[[i]]$signal[i0]
    }
    i2th[k1 + 1] <- i
    iq2th[k1 + 1] <- ij[(col(ij) == (n+i)) & (row(ij) == (n+i))]
    k1 <- k1 + 1
    i0 <- which(el[[i]]$parent)
    nj <- length(i0)
    if(nj>0) {
      j0 <- el[[i]]$id[i0]
      i1th[k2 + 1:nj] <- j0
      sth[k2 + 1:nj] <- el[[i]]$signal[i0] ## carry on the signal
      j <- n + j0
      iq1th[k2 + 1:nj] <- ij[, n+i][j]
      k2 <- k2 + nj
      i1th[k2 + 1:nj] <- j0
      sth[k2 + 1:nj] <- el[[i]]$signal[i0] ## carry on the signal
      iq1th[k2 + 1:nj] <- ij[n+i, ][j]
      k2 <- k2 + nj
    }
  }
  stopifnot(k1 == p)
  stopifnot(k2 == (2*(p-1)))
  return(list(
    n = as.integer(n),
    p = as.integer(p),
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
#' Function to build Q from a graph model
#' @param theta vector with the log of the parameters
#' @return precision matrix
#' @export
#' @examples
#' dag1 <- list(
#'      p1 ~ p2 - p3 + c1 + c2,
#'      p2 ~ -c3 + c4, p3 ~ c5 + c6)
#' Q1 <- dag_precision(dag1, theta = c(1, 1, 1))
#' cov1 <- chol2inv(chol(Q1))
#' round(100 * cov2cor(cov1[1:6, 1:6]))
#'
#' dag2 <- list(
#'      p1 ~ p2 + c1 + c2,
#'      p2 ~ -p3 -c3 + c4, p3 ~ c5 + c6)
#' Q2 <- dag_precision(dag2, theta = c(1, 1, 1))
#' cov2 <- chol2inv(chol(Q2))
#' round(100 * cov2cor(cov2[1:6, 1:6]))
dag_precision <- function(dag, theta, debug = FALSE, new = TRUE) {
  if(new) {
    q.el <- dag_precision_elements(dag)
    Q <- q.el$q
    nc <- q.el$n
    Q[q.el$iq2th] <- Q[q.el$iq2th] +
        exp(-2 * theta[q.el$i2th])
    Q[q.el$iq1th] <- -1.0 * q.el$sth * exp(-2*theta[q.el$i1th])
  } else {
    nS <- length(dag)
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
    for(k in 1:nS) {
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

dag_correlation_elements <- function(dag) {
   d.el <- dag_elements(dag)
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

dag_covariance <- function(dag, theta) {
  ij <- dag_correlation_elements(dag)
  np <- length(ij$iv)
  nc <- length(ij$iparent)
  vv <- diag(nc)
  for(i in 1:nc) {
    for (j in 1:nc) {
      vi <- sum(exp(theta[ij$iv[[ij$iparent[i]]]]))
      vj <- sum(exp(theta[ij$iv[[ij$iparent[j]]]]))
      vv[i, j] <- vv[i, j] + vi * vj
    }
  }
  return(vv)
}
