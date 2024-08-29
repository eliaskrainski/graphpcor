#' Functions to implement Kronecker (product) models
#' as methods for kronecker()
#' @export

setMethod(
  "kronecker",
  c(X="inla.cgeneric", Y = "inla.cgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    old <- X$old | Y$old
    if(length(old) == 0)
      old <- FALSE
    if(old) { ## not need to check! (because list(a=1,b=2,a=3,b=4) is ok)
      ## Check for duplicated data names
      ## TO DO: allow it
      stopifnot(all(!(names(X$f$cgeneric$data$ints[-(1:2)]) %in%
                        names(Y$f$cgeneric$data$ints[-(1:2)]))))
      stopifnot(all(!(names(X$f$cgeneric$data$doubles) %in%
                        names(Y$f$cgeneric$data$doubles))))
      stopifnot(all(!(names(X$f$cgeneric$data$characters[-(1:2)]) %in%
                        names(Y$f$cgeneric$data$characters[-(1:2)]))))
      stopifnot(all(!(names(X$f$cgeneric$data$matrices) %in%
                        names(Y$f$cgeneric$data$matrices))))
      stopifnot(all(!(names(X$f$cgeneric$data$smatrices) %in%
                        names(Y$f$cgeneric$data$smatrices))))
    }

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
    } else {
      debug <- eval(mcall$debug)
      stopifnot(is.logical(debug))
    }

    if(is.null(mcall$useINLAprecomp)) {
      useINLAprecomp = FALSE
    }
    if(is.null(mcall$libpath)) {
      libpath <- NULL
    }

    model <- ifelse(
      old,
      "inla_cgeneric_kronecker_old",
      "inla_cgeneric_kronecker")
    if (is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("corGraphs")
      } else {
        libpath <- system.file("libs", package = "corGraphs")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "corGraphs.dll")
        } else {
          libpath <- file.path(libpath, "corGraphs.so")
        }
      }
    }

    n1 <- as.integer(X$f$n)
    n2 <- as.integer(Y$f$n)
    N <- as.integer(n1 * n2)

    if(debug) {
      cat('n1:', n1, "n2:", n2, "n:", N, "")
    }

    ## index for the Q = Q1 (x) Q2 matrix
    ij1 <- cgeneric_get(
      cgeneric_model = X,
      cmd = "graph",
      optimize = TRUE
    )
    idx <- which(ij1[[1]]<=ij1[[2]])
    ij1 <- list(
      i = ij1[[1]][idx],
      j = ij1[[2]][idx]
    )
    M1 <- length(ij1[[1]])
    if(debug) {
      cat('M1:', M1, "")
    }

    ij2 <- cgeneric_get(
      cgeneric_model = Y,
      cmd = "graph",
      optimize = TRUE
    )
    idx <- which(ij2[[1]]<=ij2[[2]])
    ij2 <- list(
      i = ij2[[1]][idx],
      j = ij2[[2]][idx]
    )
    M2 <- length(ij2[[1]])
    if(debug) {
      cat('M2:', M2, "")
    }

    idx2e <- which(ij2$i < ij2$j)
    stopifnot(length(idx2e) == (M2-n2))
    M2e <- M2 + length(idx2e)

    ii <- rep(ij1$i * n2, each = M2e) +
      c(ij2$i, ij2$j[idx2e])
    jj <- rep(ij1$j * n2, each = M2e) +
      c(ij2$j, ij2$i[idx2e])
    iiord <- order(ii)
    ijs <- which(ii[iiord] <= jj[iiord])
    ije <- list(
      ii = ii[iiord[ijs]],
      jj = jj[iiord[ijs]],
      ord = iiord[ijs]
    )
    if(debug) {
      print(str(ije))
    }
    M <- length(ijs)
    if(debug) {
      cat('M:', M, "\n")
    }

    ne1 <- (M1-n1)
    ne2 <- (M2-n2)
    nnz1 <- n1 + ne1*2
    nnz2 <- n2 + ne2*2
    NNZ <- nnz1 * nnz2
    NNZu <- (NNZ - n1*n2)/2 + n1*n2

    stopifnot(M == NNZu)

    if(debug) {
      cat("ne1:", ne1, "ne2:", ne2,
          "nnz1:", nnz1, "nnz2:", nnz2,
          "nnz:", NNZ, "NNZu:", NNZu, "\n")
    }


    ## intial data
    ret <- list(
      f = list(
        model = "cgeneric",
        n = as.integer(N),
        cgeneric = list(
          model = model,
          shlib = libpath,
          n = as.integer(N),
          debug = as.integer(debug),
          data = list()
        )
      )
    )

    ## data size for each model
    ndata1 <- sapply(
      X$f$cgeneric$data,
      length)
    ndata2 <- sapply(
      Y$f$cgeneric$data,
      length)

    ### n and the data size info
    if(old) {
      ints <- list(
        n = as.integer(
          c(N, M,
            n1, M1, ndata1,
            n2, M2, ndata2
          )
        ),
        debug = debug
      )
      ret$f$cgeneric$data$ints <-
        c(
          ints,
          X$f$cgeneric$data$ints[-(1:2)],
          Y$f$cgeneric$data$ints[-(1:2)],
          list(
            idx2e = as.integer(idx2e - 1)
          )
        )
    } else {
      ints <- list(
        n = as.integer(
          c(n1, ndata1, M1,
                ndata2, M2, N, M)
        )
      )
      if(debug) {
        cat(ints$n, "\n")
      }
      ret$f$cgeneric$data$ints <-
        c(
          ints,
          ### concatenate ints from each model
          X$f$cgeneric$data$ints[-1], ## n for model 1 is already there
          Y$f$cgeneric$data$ints,
          list(
            idx2e = as.integer(idx2e - 1)
          )
        )
    }

    ret$f$cgeneric$data$doubles <-
      c(
        X$f$cgeneric$data$doubles,
        Y$f$cgeneric$data$doubles
      )

    ret$f$cgeneric$data$characters <-
      c(
        list(
          model = model,
          shlib = libpath
        ),
        X$f$cgeneric$data$characters,
        Y$f$cgeneric$data$characters
      )

    ret$f$cgeneric$data$matrices <-
      c(
        X$f$cgeneric$data$matrices,
        Y$f$cgeneric$data$matrices
      )

    ret$f$cgeneric$data$smatrices <-
      c(
        X$f$cgeneric$data$smatrices,
        Y$f$cgeneric$data$smatrices,
        list(Kgraph = c(
          N, N, M,
          ije$ii,
          ije$jj,
          as.numeric(ije$ord-1)
          )
        )
      )

    class(ret) <- "inla.cgeneric"
    class(ret$f$cgeneric) <- "inla.cgeneric"

    if(is.null(X$f$extraconstr)) {
      if(is.null(Y$f$extraconstr)) {
        if(debug)
          cat("No extraconstr!\n")
      } else {
        c2 <- Y$f$extraconstr
        ret$f$extraconstr <- list(
          A = kronecker(diag(X$f$n), c2$A),
          e = rep(c2$e, X$f$n)
        )
      }
    } else {
      c1 <- X$f$extraconstr
      if(is.null(Y$f$extraconstr)) {
        ret$f$extraconstr <- list(
          A = kronecker(c1$A, diag(Y$f$n)),
          e = rep(c1$e, each = Y$f$n)
        )
      } else {
        c2 <- m1$f$extraconstr
        ret$f$extraconstr <- list(
          A = rbind(
            kronecker(c1$A, diag(ncol(c2$A))),
            kronecker(diag(ncol(c1$A)), c2$A)
          ),
          e = c(rep(c1$e, each = ncol(c2$A)),
                rep(c2$e, ncol(c1$A)))
        )
      }
    }

    return(ret)

  }
)

setMethod(
  "kronecker",
  c(X="inla.cgeneric", Y = "inla.rgeneric"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

    mcall <- match.call()
    if(is.null(mcall$debug)) {
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
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

### follows INLA:::inla.rgeneric.define() but no assign env
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
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
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

### follows INLA:::inla.rgeneric.define() but no assign env
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
      debug <-
        max(X$f$cgeneric$debug,
            Y$f$cgeneric$debug)
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

