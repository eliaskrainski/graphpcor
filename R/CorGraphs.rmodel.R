#' R generic implementation of the graph model.
#' @param cmd the name of the object to return.
#' @param theta internal parameters.
#' @description
#' An argument 'args' should be provided as a named list with
#' \itemize{
#'  \item SP: output from the [GraphDens()] function.
#'  \item GraphPrior: function to compute the prior.
#'  \item S, lambda, Tdist: required to [GraphPrior()]
#' }
#' @return the asked object
#' @export
CorGraphs.rmodel <- function(cmd = c("graph", "Q", "mu", "initial",
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
    q <- exp(theta[-(1:iArgs$SP$NC)])
    Qmat <- numDeriv::hessian(iArgs$SP$JD[[1]], rep(0, length(iArgs$SP$STR[[1]])), SDev=q)*iArgs$SP$Pmat # precision matrix
    diag(Qmat) <- diag(Qmat) + 1e-12
      VC <- solve(Qmat) # variance-covariance
      LUB <- VC[1:iArgs$SP$NC,1:iArgs$SP$NC] # left-upper block
      Cor <- diag(diag(LUB)^(-1/2)) %*% LUB %*% diag(diag(LUB)^(-1/2)) # correlation
      SD <- diag(exp(-1/2*theta[1:iArgs$SP$NC]))
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
