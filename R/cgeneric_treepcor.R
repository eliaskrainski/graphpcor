#' Build an `cgeneric` for [treepcor()])
#' @description
#' This set the necessary data to implement the penalized
#' complexity prior for a correlation matrix considering
#' a three as proposed in
#' [Sterrantino et. al. 2025](https://arxiv.org/abs/2312.06289)
#' @param graph object of class `treepcor` for the model specification.
#' @param lambda the lambda parameter for the graph correlation prior.
#' @param sigma.prior.reference a vector with the reference values
#' to define the prior for the standard deviation parameters.
#' @param sigma.prior.probability a vector with the probability values
#' to define the prior for the standard deviation parameters.
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param libpath string, default is NULL, with the path to the shared object.
#' @details
#'  The correlation prior as in the paper depends on the lambda value.
#'  The prior for each \eqn{sigma_i} is the Penalized-complexity prior
#'  which can be defined from the following probability statement
#'  P(sigma > U) = a.
#' where "U" is a reference value and "a" is a probability.
#' The values "U" and probabilities "a" for each \eqn{sigma_i}
#' are passed in the `sigma.prior.reference` and `sigma.prior.probability`
#' arguments.
#' If a=0 then U is taken to be the fixed value of the corresponding sigma.
#' E.g. if there are three sigmas in the model and one supply
#'  sigma.prior.reference = c(1, 2, 3) and
#'  sigma.prior.probability = c(0.05, 0.0, 0.01)
#' then the sigma is fixed to 2 and not estimated.
#' @seealso [treepcor()] and [cgeneric()]
#' @return a `inla.cgeneric`, [cgeneric()] object.
#' @useDynLib graphpcor, .registration = TRUE
#' @export
cgeneric_treepcor <-
  function(graph,
           lambda,
           sigma.prior.reference,
           sigma.prior.probability,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {


    dd <- dim(graph)
    d.el <- edges(graph)[1:dd[2]]
    ich <- unlist(lapply(d.el, function(x)
      x$id[!x$parent]))
    sch <- unlist(lapply(d.el, function(x)
      x$sign[!x$parent]))
    sch <- sch[ich]
    if(debug) {
      cat(c(sch = sch), "\n")
    }
    d.elc <- etreepcor2variance(d.el[1:dd[2]])
    if(debug) {
      print(str(d.elc))
    }
    dd <- dim(graph)
    np <- dd[2]
    nv <- sapply(d.elc$iv, length)
    if(debug)
      cat("np = ", np, " and nv: ", nv, "\n")
    iiv <- rep(1:np, nv)
    jjv <- unlist(lapply(d.elc$iv, sort))
    itop <- d.elc$itop
    if(debug) {
      cat(c(iiv=iiv), "\n")
      cat(c(jjv=jjv), "\nitop:\n")
      print(itop)
    }
    nc <- nrow(itop)
    ii <- col(itop)[!upper.tri(itop)]
    jj <- row(itop)[!upper.tri(itop)]
    if(debug) {
      print(str(list(nc=nc,ii=ii,jj=jj)))
    }

    stopifnot(length(sigma.prior.reference) == nc)
    stopifnot(length(sigma.prior.probability) == nc)
    stopifnot(all(sigma.prior.probability>0.0))
    stopifnot(all(sigma.prior.probability<1.0))
    slambdas <- -log(sigma.prior.probability) / sigma.prior.reference

    stopifnot(lambda>0)

    the_model <- cgeneric.default(
      model = "inla_cgeneric_treepcor",
      libpath = libpath,
      debug = as.logical(debug),
      useINLAprecomp = as.logical(useINLAprecomp),
      n = as.integer(nc),
      np = as.integer(np),
      nv = as.integer(nv),
      ipar = as.integer(d.elc$iparent-1L),
      iiv = as.integer(iiv-1L),
      jjv = as.integer(jjv-1L),
      itop = as.integer(itop-1L),
      ii = as.integer(ii-1L),
      jj = as.integer(jj-1L),
      lambda = as.double(lambda),
      slambdas = as.double(slambdas),
      schildren = as.double(sch)
    )

    return(the_model)

  }
