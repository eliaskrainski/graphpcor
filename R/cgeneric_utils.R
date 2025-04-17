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
