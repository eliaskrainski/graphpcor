#'
#'  R generic implementation of the graph model
#'  @export
CorGraphs.model <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                    "log.prior", "quit"), theta = NULL, Argm = NULL)
{
  envir <- parent.env(environment())
  if (exists("GS", envir=envir)) {
    GS <- get("GS", envir=envir) # need to get or already in environment?
  } else {
    GS <- Argm$SP # graph structure
    assign("GS", GS, envir=envir)
  }

  graph <- function() {
    return(Q()) # does not depend on theta
  }

  Q <- function() {
    q <- exp(theta[-(1:Argm$SP$NC)])
    Q <- hessian(GS$JD[[1]], rep(0, length(GS$STR[[1]])), SDev=q)*GS$Pmat # precision matrix
    diag(Q) <- diag(Q)+1e-12
    VC <- solve(Q) # variance-covariance
    LUB <- VC[1:Argm$SP$NC,1:Argm$SP$NC] # left-upper block
    Cor <- diag(diag(LUB)^(-1/2)) %*% LUB %*% diag(diag(LUB)^(-1/2)) # correlation
    SD <- diag(exp(-1/2*theta[1:Argm$SP$NC]))
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
    q <- theta[-(1:Argm$SP$NC)]
    val <- (sum(stats::dgamma(exp(theta[1:Argm$SP$NC]), shape = 1, rate = 0.01, log = TRUE)) + sum(theta[1:Argm$SP$NC]) +
              sum(unlist(Argm$GraphPrior(Argm$S, q, Argm$lambda, Argm$SP, Argm$Tdist)))+sum(theta[-(1:Argm$SP$NC)]))
    return(val)
  }

  initial <- function() {
    return(c(rep(4, Argm$SP$NC), rep(Argm$init, Argm$SP$NP)))
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
