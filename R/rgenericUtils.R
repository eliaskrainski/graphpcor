#' Define rgeneric methods.
#' @rdname rgeneric
#' @param model an object used to define the model.
#' Its class will define which method is considered.
#' @param debug logical indicating debug state.
#' @param compile logical indicating to compile the model.
#' @param optimize logical indicating if only the elements
#' of the precision matrix are returned.
#' @param ... additional arguments to be used internally
#' for the model, for example, additional data.
#' @export
rgeneric <- function(model,
                     debug = FALSE,
                     compile = TRUE,
                     optimize = TRUE,
                     ...) {
  UseMethod("rgeneric")
}
#' The rgeneric default method.
#' @rdname rgeneric
#' @param model the model defined as a function.
#' See the 'rgeneric' vignette from the INLA package.
#' @export
#' @importFrom INLA inla.cgeneric.define
rgeneric.default <- function(model,
                             debug = FALSE,
                             compile = TRUE,
                             optimize = TRUE,
                             ...) {
## it uses INLA::inla.rgeneric.define()

  rmodel <- INLA::inla.rgeneric.define(
    model = model,
    debug = debug,
    compile = compile,
    optimize = optimize,
    ...
  )

  class(rmodel) <- "inla.rgeneric"
  class(rmodel$f$rgeneric) <- "inla.rgeneric"
  return(rmodel)
}
#' The graph method for 'inla.rgeneric'
#' @export
graph.inla.rgeneric <- function(model) {
  return(INLA::inla.rgeneric.q(
      rmodel = model,
      cmd = "graph"
  ))
}
#' The precision method for an 'inla.rgeneric' object.
#' @param ... additional parameter such as 'theta'
#' If 'theta' is not supplied, initial will be taken.
#' @export
precision.inla.rgeneric <- function(model, ...) {
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
  } else {
    warning("Using the 'default' initial parameter:")
    theta <- initial(model)
    cat(theta, '\n')
  }
  INLA::inla.rgeneric.q(
    rmodel = model,
    cmd = "Q",
    theta = theta)
}
#' The initial method for 'inla.rgeneric'
#' @export
initial.inla.rgeneric <- function(model) {
  INLA::inla.rgeneric.q(
    rmodel = model,
    cmd = "initial")
}

#' The mu method for 'inla.rgeneric'
#' @export
mu.inla.rgeneric <- function(model, theta) {
  if(missing(theta)) {
    theta <- initial(model)
    warning(paste(
      "Using the default initial theta as:",
      format(theta)))
  }
  INLA::inla.rgeneric.q(
    rmodel = model,
    cmd = "mu",
    theta = theta)
}

#' The prior metho for 'inla.rgeneric'
#' @param theta the parameter.
#' @export
prior.inla.rgeneric <- function(model, theta) {
  return(INLA::inla.rgeneric.q(
    rmodel = model,
    cmd = "log.prior",
    theta = theta))
}
