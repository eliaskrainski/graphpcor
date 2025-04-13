#' @describeIn old_corGraph
#' R generic implementation of the graph model.
#' @param cmd an string to specify which model element to get
#' @param theta numeric vector with the model parameters.
#' @description
#' An argument 'args' should be provided as a named list with
#' \itemize{
#'  \item SP: output from the [GraphDens()] function.
#'  \item GraphPrior: function to compute the prior.
#'  \item S, lambda, Tdist: required to [GraphPrior()]
#'  }
#' @importFrom stats dgamma
#' @return an `inla.rgeneric`
old_rgeneric_CorGraph <- function(cmd = c("graph", "Q", "mu", "initial",
                                    "log.norm.const", "log.prior", "quit"),
                               theta = NULL)
{

  envir <- parent.env(environment())
  if (!exists("iArgs", envir = envir)) {
    assign("iArgs", args, envir = envir)
  }
  iArgs <- get("iArgs", envir = envir)

  graph <- function() {
    return(Q()) # does not depend on theta
  }

  Q <- function() {
    if(iArgs$new) {
      S <- iArgs$S
      debug = FALSE
      nS <- length(S)
      stilde <- sapply(S, function(x)
        strsplit(as.character(x), split = "~"))
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
      elements2 <- lapply(strsplit(
        gsub(" ", "", unlist(stilde[3, ])),
        "+", fixed = TRUE), unique)
      if(debug)
        print(elements2)
      stopifnot(all(unique(
        unlist(lapply(elements2, substr, 1, 1))
      ) %in% c("p", "c")))
      iielements2 <- lapply(elements2, function(x)
        as.integer(substring(x, 2)))
      if(debug)
        print(iielements2)
      iparent2 <- lapply(elements2, function(x)
        substr(x, 1, 1) == "p")
      if(debug)
        print(iparent2)
      stopifnot(all(unlist(
        iielements2)[unlist(iparent2)]) %in% iParents1)
      ichildren2 <- lapply(elements2, function(x)
        substr(x, 1, 1) == "c")
      if(debug)
        print(ichildren2)
      NC <- sum(unlist(ichildren2))
      if(debug)
        cat("NC =", NC, "\n")
      stopifnot(max(unlist(
        iielements2)[unlist(ichildren2)]) == NC)

        q2i <- exp(theta[-(1:iArgs$SP$NC)])
        q2i

        Q <- cbind(rbind(diag(NC), matrix(0, NP, NC)),
                   rbind(matrix(0, NC, NP), diag(q2i)))
        Q
        jj <- NC + iParents1
        for(k in 1:nS) {
          if (any(ichildren2[[k]])) {
            ii <- iielements2[[k]][ichildren2[[k]]]
            cat("k1 =", k, "ii =", ii, "\n")
            Q[ii, jj[k]] <- -1
            Q[jj[k], jj[k]] <- Q[jj[k], jj[k]] + length(ii)
          }
          if (any(iparent2[[k]])) {
            i1 <- NC + iParents1[k]
            ii <- NC + iielements2[[k]][iparent2[[k]]]
            cat("k2 =", k, "ii =", ii, "\n")
            for(i in ii) {
              Q[i1, i1] <- Q[i1, i1] + q2i[i-NC]
              Q[i, i1] <- Q[i, i1] - q2i[i-NC]
            }
          }
        }
        Qmat <- Q

    } else {
    q <- exp(theta[-(1:iArgs$SP$NC)])
    Qmat <- numDeriv::hessian(iArgs$SP$JD[[1]], rep(0, length(iArgs$SP$STR[[1]])), SDev=q)*iArgs$SP$Pmat # precision matrix
    }
    diag(Qmat) <- diag(Qmat) + 1e-12
      VC <- solve(Qmat) # variance-covariance
      LUB <- VC[1:iArgs$SP$NC,1:iArgs$SP$NC] # left-upper block
      Cor <- diag(diag(LUB)^(-1/2)) %*% LUB %*% diag(diag(LUB)^(-1/2)) # correlation
      SD <- diag(exp(-0.5*theta[1:iArgs$SP$NC]))
      COV <- SD %*% Cor %*% SD
      return (solve(COV))

  }

  mu <- function() {
    return(numeric(0))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    q <- theta[-(1:iArgs$SP$NC)]
    val <- (sum(stats::dgamma(exp(theta[1:iArgs$SP$NC]), shape = 1, rate = 0.01, log = TRUE)) + sum(theta[1:iArgs$SP$NC]) +
              sum(unlist(iArgs$GraphPrior(iArgs$S, q, iArgs$lambda, iArgs$SP, iArgs$Tdist)))+sum(theta[-(1:iArgs$SP$NC)]))
    return(val)
  }

  initial <- function() {
    return(c(rep(4, iArgs$SP$NC), rep(iArgs$init, iArgs$SP$NP)))
  }

  quit <- function() {
    return(invisible())
  }

  if (!length(theta)) {
    theta <- initial()
  }

  val <- do.call(match.arg(cmd), args = list())

  return(val)

}
