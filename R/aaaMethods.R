#' precision method for inla object
#' @param x the inla output
#' @param ... additional arguments (not used yet)
#' @export
precision.inla <- function(x, ...) {
  if(is.null(x$misc$config$config)) {
    warning("Running inla(..., control.compute = list(..., config = TRUE))!")
    x$.args$control.compute$config <- TRUE
    x <- do.call("inla", args = x$.args)
  }
  Qu <- inla.as.sparse(
    x$misc$config$config[[1]]$Qprior
  )
#  ii <- which(Qu@i < Qu@j)
 # if(length(ii)>0) {
    Q <- #inla.as.sparse(
      sparseMatrix(
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
