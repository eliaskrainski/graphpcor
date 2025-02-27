#' Directed Tree Graph - DTG
#' @rdname dtg
setClass("dtg")

#' corgraph: Correlation graph mode
#' @rdname corgraph
setClass("corgraph")

#' GMRF model definition using the rgeneric in INLA
#' @rdname rgeneric
setClass(
  "inla.rgeneric",
  slots = "f",
  validity = function(object) {
    all(c("model", "n", "rgeneric") %in%
          names(object$f))
  }
)

