#' setClass

setClass(
  "inla.cgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object))
  }
)

setClass(
  "inla.rgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object))
  }
)

