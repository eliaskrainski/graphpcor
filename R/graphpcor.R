#' Build the objects needed to implement a
#' models based on graphs. 
#' Prepare objects to be used as a model 
#' in a `INLA` `f()` model component.
#' @useDynLib graphpcor, .registration = TRUE
#' @importFrom Matrix t
#' @export
graphpcor <- function() {
  cat("Welcome to the 'graphpcor' package!\n")
}
