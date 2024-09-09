#' Define cgeneric methods
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

#' @export
graph <- function(x, ...) {
  UseMethod("graph")
}
#' @export
graph.inla.cgeneric <- function(x, ...) {
  args <- list(...)
  if(any(names(args) == "optimize")) {
    return(cgeneric_get(x, "graph", ...))
  } else {
    return(cgeneric_get(x, "graph", optimize = FALSE))
  }
}

#' @export
precision.inla.cgeneric <- function(x, ...) {
  mc <- list(...)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
  } else {
    theta <- initial(x)
  }
  if(any(nargs == "optimize")) {
    optimize <- mc$optimize
  } else {
    optimize <- FALSE
  }
  stopifnot(is.logical(optimize))
  cgeneric_get(x, cmd = "Q", theta = theta, optimize = optimize)
}

#' @export
initial <- function(x) {
  UseMethod("initial")
}
#' @export
initial.inla.cgeneric <- function(x) {
  cgeneric_get(x, "initial")
}


#' @export
mu <- function(x) {
  UseMethod("mu")
}
#' @export
mu.inla.cgeneric <- function(x) {
  cgeneric_get(x, "mu")
}

#' Define prior methods.
#' @param model the model
#' @param additional arguments
#' @param theta a numeric vector with the model parameters.
#' @export
prior <- function(model, theta) {
  UseMethod("prior")
}
#' Prior for the 'inla.cgeneric' model.
#' @export
prior.inla.cgeneric <- function(model, theta) {
  return(cgeneric_get(cmodel = model,
                      cmd = "log_prior",
                      theta = theta))
}
#' Function to extract cgeneric model
#' @param cmodel an object containing the cgeneric model
#' @param optimize logical indicating if the graph and Q are
#' returned only the elements (if TRUE) or to be built (if FALSE).
#' If NULL (default) will use the initial from the cgeneric model.
#' @param cmd an string to specify which model element to get
cgeneric_get <- function(cmodel,
                         cmd = c("graph", "Q", "initial", "mu", "log_prior"),
                         theta,
                         optimize = TRUE
                         ) {

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

##  print(c(cmd = cmd))

  cgdata <- cmodel$f$cgeneric$data
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
    }
  }

  if(length(cmd) == 1) {
    ret <- .Call(
      "cgeneric_element_get",
      cmd,
      theta,
      cgdata$ints,
      cgdata$doubles,
      cgdata$characters,
      cgdata$matrices,
      cgdata$smatrices,
      PACKAGE = "corGraphs"
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
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "corGraphs"
        )
      }
      ret <- sparseMatrix(
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
          cgdata$ints,
          cgdata$doubles,
          cgdata$characters,
          cgdata$matrices,
          cgdata$smatrices,
          PACKAGE = "corGraphs"
        )
      }
    )
  if(optimize) {
     return(ret)
  }

  if(any(cmd == "graph")) {
    ret$graph <-
      sparseMatrix(
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
        cgdata$ints,
        cgdata$doubles,
        cgdata$characters,
        cgdata$matrices,
        cgdata$smatrices,
        PACKAGE = "corGraphs"
      )
      ij <- sparseMatrix(
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
