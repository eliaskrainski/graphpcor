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
initial <- function(x) {
  UseMethod("initial")
}
#' @export
initial.inla.cgeneric <- function(x) {
  cgeneric_get(x, "initial")
}

#' @export
precision.inla.cgeneric <- function(x, ...) {
  mc <- lapply(
    match.call(
      expand.dots = TRUE)[-1],
    eval)
  nargs <- names(mc)
  if(any(nargs == "theta")) {
    theta <- mc$theta
  } else {
    theta <- initial(x)
  }
  if(any(nargs == "optimize")) {
    optimize <- mc$optimize
  } else {
    optimize = TRUE
  }
  cgeneric_get(x, theta = theta, optimize = optimize)
}
#' Function to extract cgeneric model
#' @param cgeneric_model an object containing the cgeneric model
#' @param theta a numeric vector with theta.
#' @param optimize logical indicating if the graph and Q are
#' returned only the elements (if TRUE) or to be built (if FALSE).
#' If NULL (default) will use the initial from the cgeneric model.
#' @param cmd an string to specify which model element to get
cgeneric_get <- function(cgeneric_model,
                         cmd = c("graph", "Q", "initial", "mu", "log.prior"),
                         theta = NULL,
                         optimize = TRUE
                         ) {

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

##  print(c(cmd = cmd))

  cgdata <- cgeneric_model$f$cgeneric$data
  stopifnot(!is.null(cgdata))
  stopifnot(!is.null(cgdata$ints))
  stopifnot(!is.null(cgdata$characters))

  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

  if(is.null(theta) &
     any(cmd == c("Q", "log_prior"))) {
    theta <- .Call(
      "cgeneric_element_get",
      "initial",
      NULL,
      cgdata$ints,
      cgdata$doubles,
      cgdata$characters,
      cgdata$matrices,
      cgdata$smatrices,
      PACKAGE = "corGraphs"
      )
    if((length(cmd) == 1) & (cmd == "initial")) {
      return(theta)
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
    if(optimize | any(cmd == c("mu", "log_prior"))) {
      return(ret)
    }
    if(cmd == "graph") {
      ij <- ret
      ret <- rep(1, length(ij[[1]]))
      ##print(str(list(ij=ij, ret=ret, cmd = cmd)))
      ##print(sapply(list(ij=ij, ret=ret, cmd = cmd), summary))
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
      ##print(str(list(ij=ij, ret=ret, e=2)))
    }
    idx <- which(ij[[1]] <= ij[[2]])
    if(length(idx)>0) {
      ij <- list(
        ij[[1]][idx],
        ij[[2]][idx]
      )
      ret <- ret[idx]
    }
    return(
      sparseMatrix(
        i = ij[[1]] + 1L,
        j = ij[[2]] + 1L,
        x = ret,
        symmetric = TRUE,
        repr = "T"
      )
    )
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
