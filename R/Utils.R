#' Function that returns density for a given model structure
#' @param S model structure given as a formula list
#' @return a list with same length of S.
#' Each element as a named list to specify the left hand side
#' of the graph terms.
#' \itemize{
#'  \item parent: a logical vector indicating if the term is parent
#'  \item id: the number identifying the term
#'  \item signal: either minus or plus one
#'  }
#' @export
#' @examples
#'  S <- list(
#'      p3 ~ p1 - p2,
#'      p1 ~ -c1 + c2 + c3, p2 ~ c4)
#' str(Sleft(S))
Sleft <- function(S) {
  nS <- length(S)
  stilde <- sapply(S, function(x)
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

  stilde

  S.elements <- vector("list", nS)
  names(S.elements) <- stilde[2, ]
  for(k in 1:nS) {
    ch <- stilde[3, k]
    if(substr(ch, 1, 1) != "-")
      ch <- paste0("+", ch)
    p <- s <- integer(nchar(ch))
    i1 <- gregexpr("-", ch, fixed = TRUE)[[1]]
    if(any(i1>0)) {
      s[i1[i1>0]] <- -1
      p[i1[i1>0]] <- i1[i1>0]
    }
    i2 <- gregexpr("+", ch, fixed = TRUE)[[1]]
    if(any(i2>0)) {
      s[i2[i2>0]] <- 1
      p[i2[i2>0]] <- i2[i2>0]
    }
    ss <- strsplit(substring(ch, 2),
                   split = "+", fixed = TRUE)[[1]]
    ss <- unlist(lapply(ss, function(x)
      strsplit(x, "-", fixed = TRUE)[[1]]))
    S.elements[[k]] <- list(
      n = length(ss),
      term = ss,
      parent = substr(ss, 1, 1) == "p",
      id = as.integer(substring(ss, 2)),
      signal = s[s!=0]
    )
  }
  return(S.elements)
}
#' Function to build Q from a graph model
#' @param theta vector with the log of the parameters
#' @examples
#'  S <- list(
#'      p3 ~ p1 - p2,
#'      p1 ~ -c1 + c2 + c3, p2 ~ c4)
#' Q <- QS(S)
QS <- function(S, debug = FALSE) {
  nS <- length(S)
  stilde <- sapply(S, function(x)
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

  NC <- sum(sapply(S.elements, function(s)
    sum(!s$parent)))

  q <- exp(theta)
  S.elements <- Sleft(S)

  Qa <- cbind(rbind(diag(NC), matrix(0, NP, NC)),
              rbind(matrix(0, NC, NP), diag(q)))
  Qa
  jj <- NC + iParents1
  for(k in 1:nS) {
    Sk <- S.elements[[k]]
    ich <- !Sk$parent
    if (any(ich)) {
      ii <- Sk$id[ich]
      Qa[jj[k], jj[k]] <- Qa[jj[k], jj[k]] + length(ii)
      Qa[ii, jj[k]] <- -Sk$signal[ich]
    }
    if (any(Sk$parent)) {
      i1 <- NC + iParents1[k]
      ii <- NC + Sk$id[Sk$parent]
      sq <- -Sk$signal[Sk$parent] * q[ii-NC]
      for(i in 1:length(ii)) {
        Qa[i1, i1] <- Qa[i1, i1] + q[ii[i] - NC]
        Qa[ii[i], i1] <- Qa[ii[i], i1] + sq[i]
      }
    }
  }
  Qa

  Q <- Qa
  Q[lower.tri(Q)] <- t(Qa)[lower.tri(Q)]
  return(Q)
}

