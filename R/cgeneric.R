#' Define a `inla.cgeneric' object used to define a
#' latent GMRF model using the C interface for `INLA`.
#' @rdname cgeneric
#' @seealso [INLA::cgeneric()]
#' @description
#' This organize data needed on the C interface
#' for building latent models, which are characterized
#' from a given model parameters \eqn{\theta} and the
#' the following model elements.
#'  *  `graph` to define the non-zero precision matrix pattern.
#'  only the upper triangle including the diagonal is needed.
#'  The order should be by line.
#'  * `Q` vector where the
#'     * first element (N) is the size of the matrix,
#'     * second element (M) is the number of non-zero
#'     elements in the upper part (including) diagonal
#'     * the remaining (M) elements are the actual
#'     precision (upper triangle plus diagonal) elements
#'     whose order shall follow the graph definition.
#'  * `mu` the mean vector,
#'  * `initial` vector with
#'    * first element as the number of the parameters in the model
#'    * remaining elements should be the initials for the model parameters.
#'  * `log.norm.const` log of the normalizing constant.
#'  * `log.prior` log of the prior for the model parameters.
#'
#' See details in [INLA::cgeneric()]
#' @param model either a `inla.cgeneric` to be used by one of the
#' Functions (see Usage and Functions sections) or
#' an object class for what a `inla.cgeneric` method exists.
#' E.g., if it is a character, a specific function
#' will be called, see the example in Methods section.
#' @return a `inla.cgeneric`, [cgeneric()] object.
#' @export
cgeneric <- function(model, ...) {
  UseMethod("cgeneric")
}
#' @describeIn cgeneric
#' This calls [INLA::inla.cgeneric.define()]
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param libpath string, default is NULL, with the path to the shared object.
#' @param ... arguments passed on.
#' @export
cgeneric.default <- function(model,
                             debug = FALSE,
                             useINLAprecomp = TRUE,
                             libpath = NULL,
                             ...) {
  ## it uses INLA::inla.cgeneric.define()
  if(is.null(libpath)) {
    if (useINLAprecomp) {
      shlib <- INLA::inla.external.lib("graphpcor")
    } else {
      libpath <- system.file("libs", package = "graphpcor")
      if (Sys.info()["sysname"] == "Windows") {
        shlib <- file.path(libpath, "graphpcor.dll")
      } else {
        shlib <- file.path(libpath, "graphpcor.so")
      }
    }
  } else {
    shlib <- libpath
  }

  args <- list(...)
  nargs <- names(args)
  if(any(nargs == ""))
    stop("Please name the arguments!")
  cmodel <- do.call( ## TO DO: make it independent of INLA:::inla.cgeneric.define
    "inla.cgeneric.define",
    c(list(model = model,
           debug = debug,
           shlib = shlib),
      list(...))
  )
  return(cmodel)
}
#' @describeIn cgeneric
#' Method for when `model` is a character.
#' E.g. cgeneric(model = "generic0")
#' calls [cgeneric_generic0]
#' @importFrom methods existsFunction
#' @export
cgeneric.character <- function(model, ...) {
  fn <- paste0("cgeneric_", model)
  if(!existsFunction(fn)) {
    fn <- "cgeneric.default"
  }
  return(do.call(what = fn,
                 args = list(...)))
}
