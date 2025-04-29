#' @describeIn cgeneric
#' `cgeneric_get` is an internal function used by
#' `graph`, `prec`, `initial`, `mu` or `prior`
#' methods for `inla.cgeneric`
#' @param model an object `inla.cgeneric` object.
#' @param cmd an string to specify which model element to get
#' @param theta numeric vector with the model parameters.
#' If missing, the [initial()] will be used.
#' @param optimize logical. If missing or FALSE,
#' the graph and precision are as a sparse matrix.
#' If TRUE, graph only return the row/col indexes and
#' precision return only the elements as a vector.
#' @return depends on `cmd`
cgeneric_get <- function(model,
                         cmd = c("graph", "Q", "initial", "mu", "log_prior"),
                         theta,
                         optimize = TRUE) {

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

  ##  print(c(cmd = cmd))

  cgdata <- model$f$cgeneric$data
  stopifnot(!is.null(cgdata))
  stopifnot(!is.null(cgdata$ints))
  stopifnot(!is.null(cgdata$characters))

  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

  args1 <- list("cgeneric_element_get", "initial", NULL, 0L)
  tr0 <- try(do.call(".Call", c(args1, model$f$cgeneric$data)),
             silent = TRUE)
  if(inherits(tr0, 'try-error')) {

    warning('Using indirect code to extract elements!')
    if(any(cmd %in% "log_prior")) {
      warning('"log prior" is not suported without the direct code!')
    }

    ret <- cgeneric.get(model, theta)
    ret$graph <- INLA::inla.as.sparse(ret$graph)
    ret$Q <- INLA::inla.as.sparse(ret$Q)
    stopifnot(all(ret$graph@i == ret$Q@i))
    stopifnot(all(ret$graph@j == ret$Q@j))
    if(optimize) {
      o <- intersect(
        order(ret$graph@i),
        which(ret$graph@j >= ret$graph@i))
      ret$graph <- list(
        ret$graph@i[o],
        ret$graph@j[o]
      )
      ret$Q@x <- ret$Q@x[o]
    }
    names(ret) <- gsub("theta", "initial", names(ret))
    names(ret) <- gsub("log.prior", "log_prior",
                       names(ret), fixed = TRUE)
    ret <- ret[cmd]
    if(length(ret) == 1) {
      ret <- ret[[1]]
    }

  } else {
    if(missing(theta)) {
      if(cmd %in% c("Q", "log_prior")) {
        stop("Please provide 'theta'!")
      } else {
        theta <- NULL
        ntheta = 0L
      }
    } else {
      if(inherits(theta, "matrix")) {
        ntheta <- as.integer(ncol(theta))
      } else {
        ntheta <- 1L
      }
    }

    if(length(cmd) == 1) {
      ret <- .Call(
        "cgeneric_element_get",
        cmd,
        theta,
        as.integer(ntheta),
        cgdata$ints,
        cgdata$doubles,
        cgdata$characters,
        cgdata$matrices,
        cgdata$smatrices,
        PACKAGE = "graphpcor"
      )

      if((cmd %in% c("graph", "Q")) && (!optimize)) {
        if(cmd == "graph") {
          ij <- ret
          ret <- rep(1, length(ij[[1]]))
        } else {
          ij <- .Call(
            "cgeneric_element_get",
            "graph",
            NULL,
            as.integer(ntheta),
            cgdata$ints,
            cgdata$doubles,
            cgdata$characters,
            cgdata$matrices,
            cgdata$smatrices,
            PACKAGE = "graphpcor"
          )
        }
        ret <- Matrix::sparseMatrix(
          i = ij[[1]] + 1L,
          j = ij[[2]] + 1L,
          x = ret,
          symmetric = TRUE,
          repr = "T"
        )
      }
      return(ret)
    }

    names(cmd) <- cmd
    ret <-
      lapply(
        cmd, function(x) {
          .Call(
            "cgeneric_element_get",
            x,
            theta,
            as.integer(ntheta),
            cgdata$ints,
            cgdata$doubles,
            cgdata$characters,
            cgdata$matrices,
            cgdata$smatrices,
            PACKAGE = "graphpcor"
          )
        }
      )
    if(optimize) {
      return(ret)
    }

    if(any(cmd == "graph")) {
      ret$graph <-
        Matrix::sparseMatrix(
          i = ret$graph[[1]] + 1L,
          j = ret$graph[[2]] + 1L,
          x = rep(1, length(ret$graph[[1]])),
          symmetric = TRUE,
          repr = "T"
        )
    }

    if(any(cmd == "Q")) {
      if(any(cmd == "graph")) {
        ij <- ret$graph
      } else {
        ij <- .Call(
          "cgeneric_element_get",
          "graph",
          theta,
          as.integer(ntheta),
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "graphpcor"
        )
        ij <- Matrix::sparseMatrix(
          i = ij[[1]] + 1L,
          j = ij[[2]] + 1L,
          symmetric = TRUE,
          repr = "T"
        )
      }
      ij@x <- ret$Q
      ret$Q <- ij
    }

  }

  return(ret)

}
#' @describeIn cgeneric
#' Retrieve the initial model parameter(s).
#' @export
initial <- function(model) {
  UseMethod("initial")
}
#' @describeIn cgeneric
#' Retrive the initial parameter(s) of an `inla.cgeneric` model.
#' @export
initial.inla.cgeneric <- function(model) {
  cgeneric_get(model, "initial")
}
#' @describeIn cgeneric
#' Evaluate the mean.
#' @export
mu <- function(model, theta) {
  UseMethod("mu")
}
#' @describeIn cgeneric
#' Evaluate the mean for an `inla.cgeneric` model.
#' @export
mu.inla.cgeneric <- function(model, theta) {
  cgeneric_get(model, "mu", theta = theta)
}
#' @describeIn cgeneric
#' Evaluate the log-prior.
#' @export
prior <- function(model, theta) {
  UseMethod("prior")
}
#' @describeIn cgeneric
#' Evaluate the prior for an `inla.cgeneric` model
#' @export
prior.inla.cgeneric <- function(model, theta) {
  return(cgeneric_get(model = model,
                      cmd = "log_prior",
                      theta = theta))
}
#' @describeIn cgeneric
#' Retrieve the graph
#' @export
graph <- function(model, ...) {
  UseMethod("graph")
}
#' @describeIn cgeneric
#' Retrieve the graph of an `inla.cgeneric` object
#' @export
graph.inla.cgeneric <- function(model, ...) {
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "optimize")) {
    optimize <- mc$optimize
  } else {
    optimize <- FALSE
  }
  stopifnot(is.logical(optimize))
  return(cgeneric_get(
    model, "graph",
    optimize = optimize))
}
#' @describeIn cgeneric
#' Evaluate [prec()] on a model
Q <- function(model, ...) {
  UseMethod("prec")
}
#' @describeIn cgeneric
#' Evaluate [prec()] on an `inla.cgeneric` object
#' @export
prec.inla.cgeneric <- function(model, ...) {
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
  } else {
    warning("Using the 'default' initial parameter:")
    theta <- initial(model)
    cat(theta, '\n')
  }
  if(any(nargs == "optimize")) {
    optimize <- mc$optimize
  } else {
    optimize <- FALSE
  }
  stopifnot(is.logical(optimize))
  cgeneric_get(model, cmd = "Q", theta = theta, optimize = optimize)
}
#' @describeIn prec R level function to
#' extract elements calling INLA program
cgeneric.get <- function(model, theta) {
  result <- INLA::inla.cgeneric.q(model)
  if(!missing(theta)) {
    ctrlf <- list(hyper = list(prec = list(
      initial = 10, fixed = TRUE
    )))
    tr0 <- try(INLA::inla(
        formula = y ~ -1 + f(one, model = model),
        data = data.frame(y = NA, one = 1),
        control.family = ctrlf,
        control.mode = list(
          restart = FALSE,
          fixed = TRUE
        ),
        silent = 2L, verbose = FALSE),
        silent = TRUE)
    if(inherits(tr0, "try-error"))
      return(tr0)
    stopifnot(length(theta) == length(tr0$mode$theta))
    iout <- INLA::inla(
      formula = y ~ -1 + f(one, model = model),
      data = data.frame(y = NA, one = 1),
      control.family = ctrlf,
      control.mode = list(
        theta = theta,
        restart = FALSE,
        fixed = TRUE
      ),
      control.compute = list(config = TRUE),
      control.inla = list(int.strategy = "eb"),
      silent = 2L, verbose = FALSE)
    result$Q <- prec(iout)
  }
  return(result)
}
