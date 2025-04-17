#' Build an `inla.cgeneric` to implement the PC-prior of a
#' precision matrix as inverse of a correlation matrix.
#' @param n integer to define the size of the matrix
#' @param lambda numeric (positive), the penalization rate parameter
#' @param theta.base numeric vector with the model parameters
#' at the base model
#' @param debug integer, default is zero, indicating the verbose level.
#' Will be used as logical by INLA.
#' @param useINLAprecomp logical, default is TRUE, indicating if it is to
#' be used the shared object pre-compiled by INLA.
#' This is not considered if 'libpath' is provided.
#' @param libpath string, default is NULL, with the path to the shared object.
#' @details
#' The precision matrix parametrization
#' step 1:
#' \deqn{Q0 = \left[
#' \begin{array}{ccccc}
#'   1 & & & & \\
#'   \theta_1 & 1 & & & \\
#'   \theta_2 & \theta_n & & & \\
#'   \vdots & & \ldots & \ddots & \\
#'   \theta_{n-1} & \theta_{2n-3} \ldots & \theta_m & 1
#'   \end{array}
#'   \right] }
#'
#' step 2: \eqn{V = Q0^{-1}}
#'
#' step 3: \eqn{S = diag(V)^{1/2}}
#'
#' step 4: \eqn{C = SVS}
#'
#' step 5: \eqn{Q = C^{-1}}
#'
#' \deqn{p(Q|\lambda) = p(\theta[1:m] | lambda) =}
#' \deqn{    p_C(C(Q)) | Jacobian C(Q) |}
#'  where p_C is the PC-prior for correlation,
#'   see section 6.2 of Simpson et. al. (2017),
#' which is based on the hypersphere decomposition.
#'
#' The hypershere decomposition, as proposed in
#' Rapisarda, Brigo and Mercurio (2007)
#' consider \eqn{\theta[k] \in [0, \infty], k=1,...,m=n(n-1)/2}
#' compute \eqn{x[k] = pi/(1+exp(-theta[k]))}
#' organize it as a lower triangle of a \eqn{n \times n} matrix
#' \deqn{B[i,j] = \left\{\begin{array}{cc}
#' cos(x[i,j]) & j=1 \\
#' cos(x[i,j])prod_{k=1}^{j-1}sin(x[i,k]) &  2 <= j <= i-1 \\
#' prod_{k=1}^{j-1}sin(x[i,k])  & j=i \\
#' 0 &  j+1 <= j <= n \end{array}\right.}
#' Result
#' \deqn{\gamma[i,j] = -log(sin(x[i,j]))}
#'  \deqn{KLD(R) = \sqrt(2\sum_{i=2}^n\sum_{j=1}^{i-1} \gamma[i,j]}
#' @references
#' Daniel Simpson, H\\aa vard Rue, Andrea Riebler, Thiago G.
#' Martins and Sigrunn H. S\\o rbye (2017).
#' Penalising Model Component Complexity:
#' A Principled, Practical Approach to Constructing Priors
#' Statistical Science 2017, Vol. 32, No. 1, 1–28.
#' <doi 10.1214/16-STS576>
#'
#' Rapisarda, Brigo and Mercurio (2007).
#'   Parameterizing correlations: a geometric interpretation.
#'   IMA Journal of Management Mathematics (2007) 18, 55-73.
#'   <doi 10.1093/imaman/dpl010>
#'
#' @return a `inla.cgeneric`, [cgeneric()] object.
cgeneric_pc_prec_correl <-
  function(n,
           lambda,
           theta.base,
           debug = FALSE,
           useINLAprecomp = TRUE,
           libpath = NULL) {

    if(is.null(libpath)) {
      if (useINLAprecomp) {
        libpath <- INLA::inla.external.lib("graphpcor")
      } else {
        libpath <- system.file("libs", package = "graphpcor")
        if (Sys.info()["sysname"] == "Windows") {
          libpath <- file.path(libpath, "graphpcor.dll")
        } else {
          libpath <- file.path(libpath, "graphpcor.so")
        }
      }
    }

    stopifnot(n>1)
    stopifnot(lambda>0)

    m <- n*(n-1)/2

    if(missing(theta.base)) {
      theta.base <- rep(0, m)
      warning("Missing 'theta.base' model. Assuming 'iid' by using:\n",
              paste(theta.base, collapse = ", "))
    }
    H.el <- theta2H(theta.base)
    if(debug) {
      cat("H elements\n")
      print(str(H.el))
    }

    ## constant: log( \lambda \pi^{m-1}/2 |H| )
    lc <- log(lambda) -(m-1)*log(pi) - log(2)
    lc <- lc - sum(log(H.el$svd$d))

    if(debug) {
      cat('log C', lc, '\n')
    }

    cmodel = "inla_cgeneric_pc_prec_correl"

    the_model <- list(
      f = list(
        model = "cgeneric",
        n = as.integer(n),
        cgeneric = list(
          model = cmodel,
          shlib = libpath,
          n = as.integer(n),
          debug = as.logical(debug),
          data = list(
            ints = list(
              n = as.integer(n),
              debug = as.integer(debug)
            ),
            doubles = list(
              lambda = as.numeric(lambda),
              lconst = as.numeric(lc)
            ),
            characters = list(
              model = cmodel,
              shlib = libpath
            ),
            matrices = list(
              ),
            smatrices = list(
              )
            )
          )
        )
      )

    class(the_model) <- "inla.cgeneric"
    class(the_model$f$cgeneric) <- "inla.cgeneric"

    return(the_model)

}
