#' Function to extract cgeneric model
#' @param cgeneric_model an object containing the cgeneric model
#' @param theta a numeric vector with theta.
#' @param optimize logical indicating if the graph and Q are
#' returned only the elements (if TRUE) or to be built (if FALSE).
#' If NULL (default) will use the initial from the cgeneric model.
#' @param cmd an string to specify which model element to get
#' @export
cgeneric_get <- function(cgeneric_model,
                         cmd = c("graph", "Q", "initial", "mu", "log.prior"),
                         theta = NULL,
                         optimize = TRUE
                         ) {

  ret <- NULL
  cmd[cmd == "log.prior"] <- "log_prior"
  cmd <- unique(cmd)

  cgdata <- cgeneric_model$f$cgeneric$data
  stopifnot(!is.null(cgdata))
  stopifnot(!is.null(cgdata$ints))
  stopifnot(!is.null(cgdata$characters))

  cmds <- c("graph", "Q", "initial", "mu", "log_prior")
  cmd <- match.arg(cmd,
                   cmds,
                   several.ok = TRUE)
  stopifnot(length(cmd)>0)

  itheta <- .Call(
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

  if(length(cmd) == 1) {
    if(cmd == "initial") {
      return(itheta)
    }
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
    if(cmd != "graph") {
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
    } else {
      ij <- ret
      ret <- rep(1, length(ij[[1]]))
    }
    return(
      INLA::inla.as.sparse(
        sparseMatrix(
          i = ij[[1]] + 1,
          j = ij[[2]] + 1,
          x = ret,
          symmetric = TRUE
        )
      )
    )
  }

  names(cmd) <- cmd
  if(is.null(theta)) {
    theta = itheta
  }
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
    }
    ret$Q <- INLA::inla.as.sparse(
      sparseMatrix(
        i = ij[[1]] + 1L,
        i = ij[[2]] + 1L,
        x = ret$Q,
        symmetric = TRUE
      )
    )
  }
  if(any(cmd == "graph")) {
    ret$graph <- INLA::inla.as.sparse(
      sparseMatrix(
        i = ret$graph[[1]] + 1L,
        i = ret$graph[[2]] + 1L,
        x = 1
      )
    )
  }
  return(ret)

}
