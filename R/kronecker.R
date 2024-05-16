#' Functions to implement Kronecker (product) models
#' as methods for kronecker()
#' @export
setMethod(
  "kronecker",
  c(X="inla.cgeneric", Y = "inla.rgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <- FALSE
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    n <- X$f$n * Y$f$n

    nth1 <- length(cgeneric_get(X, cmd = "initial"))
    nth2 <- length(inla.rgeneric.q(Y, cmd = "initial"))

    kmodel <- function(cmd = c("graph", "Q", "mu",
                               "initial", "log.norm.const",
                               "log.prior", "quit"),
                       theta = NULL) {

      graph <- function(n, theta) {
        g1 <- cgeneric_get(X, cmd = "graph")
        g1 <- sparseMatrix(
          i = g1[[1]] + 1,
          j = g1[[2]] + 1,
          x = 1,
          symmetric = TRUE
        )
        g2 <- inla.rgeneric.q(Y, cmd = "graph")
        return(kronecker(g1, g2))
      }

      Q <- function(n, theta) {
        g1 <- cgeneric_get(X, cmd = "graph")
        q1 <- cgeneric_get(X, cmd = "Q",
                           theta = theta[1:nth1])
        Q1 <- sparseMatrix(
          i = g1[[1]] + 1,
          j = g1[[2]] + 1,
          x = q1,
          symmetric = TRUE
        )
        Q2 <- inla.rgeneric.q(Y, cmd = "Q",
                              theta = theta[nth1+1:nth2])
        QQ <- inla.as.sparse(kronecker(Q1, Q2))
        idx <- which(QQ@i <= QQ@j)
        return(QQ@x[idx])
      }

      mu <- function(n, theta)
        return(numeric(0))

      log.norm.const <- function(n, theta)
        return(numeric(0))

      log.prior <- function(n, theta) {
        return(
          cgeneric_get(X, cmd = "log.prior",
                          theta = theta[1:nth1]) +
            inla.rgeneric.q(Y, cmd = "log.prior",
                         theta = theta[nth1+1:nth2])
        )
      }

      initial <- function(n, theta) {
        return(
          c(
            cgeneric_get(X, cmd = "initial"),
            inla.rgeneric.q(Y, cmd = "initial")
          )
        )
      }

      quit <- function(n, theta) {
        return(invisible())
      }

      cmd <- match.arg(cmd)

      ret <- do.call(
        cmd,
        args = list(n = n,
                    theta = theta
        )
      )

      return(ret)

    }

    ### this follows INLA:::inla.rgeneric.define() but no assign env
    rmodel <- list(
      f = list(
        model = "rgeneric",
        n = n,
        rgeneric = list(
          definition =
            compiler::cmpfun(
              kmodel,
              options = list(optimize = 3L)),
          debug = debug,
          optimize = TRUE
        )
      )
    )
    class(rmodel) <- "inla.rgeneric"
    class(rmodel$f$rgeneric) <- "inla.rgeneric"
    return(rmodel)

  }
)

setMethod(
  "kronecker",
  c(X="inla.rgeneric", Y = "inla.cgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <- FALSE
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    n <- X$f$n * Y$f$n

    nth1 <- length(inla.rgeneric.q(X, cmd = "initial"))
    nth2 <- length(cgeneric_get(Y, cmd = "initial"))

    kmodel <- function(cmd = c("graph", "Q", "mu",
                               "initial", "log.norm.const",
                               "log.prior", "quit"),
                       theta = NULL) {

      graph <- function(n, theta) {
        g1 <- inla.rgeneric.q(X, cmd = "graph")
        g2 <- cgeneric_get(Y, cmd = "graph")
        g2 <- sparseMatrix(
          i = g2[[1]] + 1,
          j = g2[[2]] + 1,
          x = 1,
          symmetric = TRUE
          )
        return(kronecker(g1, g2))
      }

      Q <- function(n, theta) {
        Q1 <- inla.rgeneric.q(X, cmd = "Q",
                              theta = theta[1:nth1])
        g2 <- cgeneric_get(Y, cmd = "graph")
        q2 <- cgeneric_get(Y, cmd = "Q",
                           theta = theta[nth1+1:nth2])
        Q2 <- sparseMatrix(
          i = g2[[1]] + 1,
          j = g2[[2]] + 1,
          x = q2,
          symmetric = TRUE
        )
        QQ <- inla.as.sparse(kronecker(Q1, Q2))
        idx <- which(QQ@i <= QQ@j)
        return(QQ@x[idx])
      }

      mu <- function(n, theta)
        return(numeric(0))

      log.norm.const <- function(n, theta)
        return(numeric(0))

      log.prior <- function(n, theta) {
        return(
          inla.rgeneric.q(X, cmd = "log.prior",
                          theta = theta[1:nth1]) +
            cgeneric_get(Y, cmd = "log.prior",
                         theta = theta[nth1+1:nth2])
        )
      }

      initial <- function(n, theta) {
        return(
          c(
            inla.rgeneric.q(X, cmd = "initial"),
            cgeneric_get(Y, cmd = "initial")
          )
        )
      }

      quit <- function(n, theta) {
        return(invisible())
      }

      cmd <- match.arg(cmd)

      ret <- do.call(
        cmd,
        args = list(n = n,
                    theta = theta
        )
      )

      return(ret)

    }

    ### this follows INLA:::inla.rgeneric.define() but no assign env
    rmodel <- list(
      f = list(
        model = "rgeneric",
        n = n,
        rgeneric = list(
          definition =
            compiler::cmpfun(
              kmodel,
              options = list(optimize = 3L)),
          debug = debug,
          optimize = TRUE
        )
      )
    )
    class(rmodel) <- "inla.rgeneric"
    class(rmodel$f$rgeneric) <- "inla.rgeneric"
    return(rmodel)

  }
)

setMethod(
  "kronecker",
  c(X="inla.rgeneric", Y = "inla.rgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <- FALSE
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    n <- X$f$n * Y$f$n

    nth1 <- length(inla.rgeneric.q(X, cmd = "initial"))
    nth2 <- length(inla.rgeneric.q(Y, cmd = "initial"))

    kmodel <- function(cmd = c("graph", "Q", "mu",
                               "initial", "log.norm.const",
                               "log.prior", "quit"),
                       theta = NULL) {

      graph <- function(n, theta) {
        g1 <- inla.rgeneric.q(X, "graph")
        g2 <- inla.rgeneric.q(Y, "graph")
        return(kronecker(g1, g2))
      }

      Q <- function(n, theta) {
        Q1 <- INLA:::inla.rgeneric.q(
          rmodel = X,
          cmd = "Q",
          theta = theta[1:nth1])
        Q2 <- INLA:::inla.rgeneric.q(
          rmodel = Y,
          cmd = "Q",
          theta = theta[nth1+1:nth2])
        QQ <- INLA::inla.as.sparse(
          kronecker(Q1, Q2))
        idx <- which(QQ@i <= QQ@j)
        return(QQ@x[idx])
      }

      mu <- function(n, theta)
        return(numeric(0))

      log.norm.const <- function(n, theta)
        return(numeric(0))

      log.prior <- function(n, theta) {
        lp1 <- inla.rgeneric.q(
          rmodel = X,
          cmd = "log.prior",
          theta = theta[1:nth1])
        lp2 <- inla.rgeneric.q(
          rmodel = Y,
          cmd = "log.prior",
          theta = theta[nth1+1:nth2])
        return(lp1 + lp2)
      }

      initial <- function(n, theta) {
        ini1 <- inla.rgeneric.q(
          rmodel = X, ## model1,
          cmd = "initial",
          theta = theta[1:nth1])
        ini2 <- inla.rgeneric.q(
          rmodel = Y, ## model2,
          cmd = "initial",
          theta = theta[nth1+1:nth2])
        return(c(ini1, ini2))
      }

      quit <- function(n, theta) {
        return(invisible())
      }

      cmd <- match.arg(cmd)

      ret <- do.call(
        cmd,
        args = list(
          n = n,
          theta = theta
        )
      )
      return(ret)
    }

    return(INLA::inla.rgeneric.define(
      model = kmodel,
      optimize = TRUE
    ))

  }
)

