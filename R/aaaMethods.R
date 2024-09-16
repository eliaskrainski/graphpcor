#' Define the variance method
#' @export
variance <- function(x, ...) {
  UseMethod("variance")
}
#' The variance default definition
#' @export
variance.default <- function(x, ...) {
  return(var(x))
}

#' Define the precision method
#' @export
precision <- function(x, ...) {
  UseMethod("precision")
}
#' The precision default method
#' @export
precision.default <- function(x, ...) {
  v <- variance(x, ...)
  return(
    forwardsolve(
      backsolve(
        chol(v)
      )
    )
  )
}

#' Define the cgeneric method
#' @param model an object used to define the model.
#' Its class will define which method is considered.
#' @param debug logical indicating debug state.
#' @param useINLAprecomp logical indicating if
#' it is to be used the shared object within INLA.
#' @param ... additional arguments to be treated
#' according to each method.
#' @export
cgeneric <- function(model, ...) {
  UseMethod("cgeneric")
}
#' @export
#' @importFrom INLA inla.cgeneric.define
cgeneric.default <- function(model,
                             debug = FALSE,
                             useINLAprecomp = TRUE,
                             ...) {
  ## it uses INLA::inla.cgeneric.define()
  if (useINLAprecomp) {
    shlib <- INLA::inla.external.lib("corGraphs")
  } else {
    libpath <- system.file("libs", package = "corGraphs")
    if (Sys.info()["sysname"] == "Windows") {
      shlib <- file.path(libpath, "corGraphs.dll")
    } else {
      shlib <- file.path(libpath, "corGraphs.so")
    }
  }

  args <- list(...)
  nargs <- names(args)
  if(any(nargs == ""))
    stop("Please name the arguments!")
  stopifnot(any(nargs == "n"))
  cmodel <- do.call(
    "inla.cgeneric.define",
    c(list(model = model,
           debug = debug,
           shlib = shlib),
      list(...))
  )
  return(cmodel)
}
#' Define the initial model to apply for a model object
#' @param model the model object
#' @export
initial <- function(model) {
  UseMethod("initial")
}

#' @export
mu <- function(model, theta) {
  UseMethod("mu")
}
#' Define prior methods.
#' @param theta a numeric vector with the model parameters.
#' @export
prior <- function(model, theta) {
  UseMethod("prior")
}
#' Define the graph generic method
#' @param ... additional arguments for each method
#' @export
graph <- function(model, ...) {
  UseMethod("graph")
}
#' Define the precision method for an inla output object
#' @param model the fitted model as an inla output
#' @export
precision.inla <- function(model, ...) {
  if(is.null(model$misc$config$config)) {
    warning("Running inla(..., control.compute = list(..., config = TRUE))!")
    model$.args$control.compute$config <- TRUE
    model <- do.call("inla", args = model$.args)
  }
  Qu <- INLA::inla.as.sparse(
    model$misc$config$config[[1]]$Qprior
  )
#  ii <- which(Qu@i < Qu@j)
 # if(length(ii)>0) {
    Q <- #inla.as.sparse(
      Matrix::sparseMatrix(
#        i = c(Qu@i, Qu@j[ii]) + 1L,
 #       j = c(Qu@j, Qu@i[ii]) + 1L,
  #      x = c(Qu@x, Qu@x[ii])
        i = Qu@i + 1L,
        j = Qu@j + 1L,
        x = Qu@x,
        symmetric = TRUE,
        repr = "T"
      )
#    )
 # } else {
  #  Q <- Qu
  #}
  return(Q)
}
