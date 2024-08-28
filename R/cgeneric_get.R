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
    cat('here\n')
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
    print(c(theta = theta))
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
      INLA::inla.as.sparse(
        sparseMatrix(
          i = ij[[1]] + 1L,
          j = ij[[2]] + 1L,
          x = ret,
          symmetric = TRUE
        )
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
        j = ij[[2]] + 1L,
        x = ret$Q,
        symmetric = TRUE
      )
    )
  }
  if(any(cmd == "graph")) {
    ret$graph <- INLA::inla.as.sparse(
      sparseMatrix(
        i = ret$graph[[1]] + 1L,
        j = ret$graph[[2]] + 1L,
        x = 1
      )
    )
  }
  return(ret)

}
