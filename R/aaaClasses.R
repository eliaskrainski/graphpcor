#' Directed Tree Graph - DTG
#' @rdname dtg
setClass("dtg")

#' Set the graph class 
#' @rdname graph
setClass("graph")

#' `inla.rgeneric` class, short `rgeneric`,
#' to define a [INLA::rgeneric()] latent model
#' @rdname rgeneric
setClass(
  "inla.rgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)

#' `inla.cgeneric` class, short `cgeneric`,
#' to define a [INLA::cgeneric()] latent model
#' @rdname cgeneric
setClass(
  "inla.cgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)

