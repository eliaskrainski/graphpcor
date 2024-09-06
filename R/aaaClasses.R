#' Directed Tree Graph - DTG.
setClass("dtg")

#' GMRF model definition using the cgeneric in INLA
setClass(
  "cgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)

#' GMRF model definition using the rgeneric in INLA
setClass(
  "rgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)

