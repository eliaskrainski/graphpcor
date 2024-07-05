#' setClass
#'

setClass(
  "extraconstr",
  slots = c("A", "e"),
  validity = function(object) {
    all(c("A", "e", "n") %in% names(object)) &&
      (ncol(object$A) == length(object$e))
  }
)

setClass(
  "inla.cgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)

setClass(
  "inla.rgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "cgeneric") %in%
          names(object$f))
  }
)

