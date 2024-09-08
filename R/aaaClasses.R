#' Directed Tree Graph - DTG.
setClass("dtg")

#' GMRF model definition using the cgeneric in INLA
setClass(
  "inla.cgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)

#' GMRF model definition using the rgeneric in INLA
setClass(
  "inla.rgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)

