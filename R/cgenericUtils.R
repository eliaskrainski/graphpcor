#' The graph method for cgeneric
#' @param model the inla.cgeneric model
#' @export
graph.inla.cgeneric <- function(model, ...) {
  args <- list(...)
  if(any(names(args) == "optimize")) {
    return(cgeneric_get(model, "graph", ...))
  } else {
    return(cgeneric_get(model, "graph", optimize = FALSE))
  }
}

#' @export
precision.inla.cgeneric <- function(model, ...) {
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

#' @export
initial.inla.cgeneric <- function(model) {
  cgeneric_get(model, "initial")
}

#' @export
mu.inla.cgeneric <- function(model) {
  cgeneric_get(model, "mu")
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
        cgdata$ints,
        cgdata$doubles,
        cgdata$characters,
        cgdata$matrices,
        cgdata$smatrices,
        PACKAGE = "corGraphs"
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
