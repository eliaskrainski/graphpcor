#' @rdname pcor-class
#' @description
#' Organize information to buidl a correlation prior
#' that penalizes the divergence from a base correlation matrix,
#' and helper functions to work with correlation matrices.
#' @param base numeric/matrix used to define the base
#' correlation matrix. If numeric vector with length 'm',
#' 'm' should be 'p(p-1)/2' in the dense case and
#' 'length(itheta)' in the sparse case.
#' @param p integer with the dimention, number of rows
#' and columns of the correlation matrix.
#' @param parametrization character to specify the
#' parametrization choice: "ITP" (or "itp"),
#' "SAP" (or "sap"), "CPC" (or "cpc"), see details.
#' @param itheta integer vector to specify the position
#' 'theta' will be placed in the (initial, before fill-in)
#' Cholesky factor in the IT. Default is missing.
#' @param d0 numeric vector to specify the diagonal
#' for the (initial) Cholesky factor in the ITP.
#' Default consider 'd0=p:1'.
#' @param decomposition character to specify which
#' decomposition is to be applied on H in order to
#' compute \eqn{\mathbf{H}^{1/2}} and
#' \eqn{\mathbf{H}^{-1/2}}. The options are
#' 'eigen', 'svd' and 'chol'
#' @param ... used to pass the arguments
#' for [numDeriv::hessian()]
#' @details
#' For 'parametrization' = "CPC" or 'parametrization' = "cpc":
#' The Canonical Partial Correlation - CPC parametrization,
#'  Lewandowski, Kurowicka, and Joe (2009), compute
#' \eqn{r[k]} = tanh(\eqn{\theta[k]}), for \eqn{k=1,...,m},
#' and the two \eqn{p\times p} matrices
#' \deqn{A = \left[
#' \begin{array}{ccccc}
#'   1 & & & & \\
#'   r_1 & 1 & & & \\
#'   r_2 & r_p & 1 & & \\
#'   \vdots & \vdots & \ddots & \ddots & \\
#'   r_{p-1} & r_{2p-3} & \ldots & r_m & 1
#' \end{array} \right]
#' \textrm{ and } B = \left[
#' \begin{array}{ccccc}
#'   1 & & & & \\
#'   \sqrt{1-r_1^2} & 1 & & & \\
#'   \sqrt{1-r_2^2} & \sqrt{1-r_p^2} & 1 & & \\
#'   \vdots & \vdots & \ddots & \ddots & \\
#'   \sqrt{1-r_{p-1}^2} & \sqrt{r_{2p-3}^2} & \ldots & \sqrt{1-r_m^2} & 1
#' \end{array} \right] }
#'
#' The matrices \eqn{A} and \eqn{B} are then used
#' to build the Cholesky factor of the correlation matrix,
#' given as
#' \deqn{L = \left[
#' \begin{array}{ccccc}
#'   1 & 0 & 0 & \ldots & 0\\
#'   A_{2,1} & B_{2,1} & 0 & \ldots & 0\\
#'   A_{3,1} & A_{3,2}B_{3,1} & B_{3,1}B_{3,2} & & \vdots \\
#'   \vdots & \vdots & \ddots & \ddots & 0\\
#'   A_{p,1} & A_{p,2}B_{p,1} & \ldots &
#'   A_{p,p-1}\prod_{k=1}^{p-1}B_{p,k} & \prod_{k=1}^{p-1}B_{p,k}
#' \end{array} \right]}
#' Note: The determinant of the correlation matriz is
#' \deqn{\prod_{i=2}^p\prod_{j=1}^{i-1}B_{i,j} = \prod_{i=2}^pL_{i,i}}

#' For 'parametrization' = "SAP" or 'parametrization' = "sap":
#' The Standard Angles Parametrization - SAP, as described in
#' Rapisarda, Brigo and Mercurio (2007), compute
#' \eqn{x[k] = \pi/(1+\exp(-\theta[k]))}, for \eqn{k=1,...,m},
#' and the two \eqn{p\times p} matrices
#' \deqn{A = \left[
#' \begin{array}{ccccc}
#'   1 & & & & \\
#'   \cos(x_1) & 1 & & & \\
#'   \cos(x_2) & \cos(x_p) & 1 & & \\
#'   \vdots & \vdots & \ddots & \ddots & \\
#'   \cos(x_{p-1}) & \cos(x_{2p-3}) & \ldots & \cos(x_m) & 1
#' \end{array} \right] \textrm{ and } B = \left[
#' \begin{array}{ccccc}
#'   1 & & & & \\
#'   \sin(x_1) & 1 & & & \\
#'   \sin(x_2) & \sin(x_p) & 1 & & \\
#'   \vdots & \vdots & \ddots & \ddots & \\
#'   \sin(x_{p-1}) & \sin(x_{2p-3}) & \ldots & \sin(x_m) & 1
#' \end{array} \right]}
#'
#' For 'parametrization' = "itp" or parametrization' = "ITP":
#' The Inverse Transform Parametrization - ITP, is applied
#' by starting with a
#' \deqn{\mathbf{L}^{(0)} = \left[ \begin{array}{ccccc}
#' \textrm{p} & 0 & 0 & \ldots & 0 \\
#' \theta_1 & \textrm{p-}1 & 0 & \ldots & 0 \\
#' \theta_2 & \theta_p & \textrm{p-}2 & \ddots & \vdots \\
#' \vdots & \vdots & \ddots & \ddots & 0 \\
#' \theta_{p-1} & \theta_{2p-3} & \ldots & \theta_m & 1
#' \end{array}\right] .}
#'
#' Then compute \eqn{\mathbf{Q}^{(0)}} =
#' \eqn{\mathbf{L}\mathbf{L}^{T}},
#' \eqn{\mathbf{V}^{(0)}} = \eqn{(\mathbf{Q}^{(0)})^{-1}} and
#' \eqn{s_{i}^{(0)}} = \eqn{\sqrt{\mathbf{V}_{i,i}^{(0)}}}, and
#' define \eqn{\mathbf{S}^{(0)}} = diag\eqn{(s_1^{(0)},\ldots,s_p^{(0)})}
#'  in order to have \eqn{\mathbf{C}}=
#'  \eqn{\mathbf{S}^{-1}\mathbf{V}^{(0)}\mathbf{S}^{-1}}.
#'
#' @references
#' Rapisarda, Brigo and Mercurio (2007).
#'   Parameterizing correlations: a geometric interpretation.
#'   IMA Journal of Management Mathematics (2007) 18, 55-73.
#'   <doi 10.1093/imaman/dpl010>
#'
#' Lewandowski, Kurowicka and Joe (2009)
#' Generating Random Correlation Matrices Based
#' on Vines and Extended Onion Method.
#' Journal of Multivariate Analysis 100: 1989–2001.
#' <doi: 10.1016/j.jmva.2009.04.008>
#'
#' Simpson, et. al. (2017)
#' Penalising Model Component Complexity:
#' A Principled, Practical Approach to Constructing Priors.
#'  Statist. Sci. 32(1): 1-28 (February 2017).
#'  <doi: 10.1214/16-STS576>
#' @returns a pcor object
#' @export
pcor <- function(base, p,
                 parametrization = 'cpc',
                 itheta, d0,
                 decomposition,
                 ...) {
  UseMethod("pcor")
}
#' @rdname pcor-class
#' @returns a `pcor` object
#' @export
#' @examples
#'
#' ## Base correlation matrix
#' c0 <- matrix(c(1,.9,.5, .9,1,.2, .5,.2,1), 3)
#'
#' ## build the 'pcor'
#' pc.c0 <- pcor(base = c0) ## base as matrix
#' pc.c0
#'
#' ## using 'theta' instead (numerically the same)
#' th0 <- pc.c0$theta
#' pc.th0 <- pcor(base = th0) ## base as vector
#' pc.th0
#'
#' all.equal(pc.c0, pc.th0)
#'
#' ## other way around
#' all.equal(pc.th0, pcor(base = pc.th0$base))
#'
#' ## ITP
#' th <- c(0.5,-1,0.5,-0.3)
#' ith <- c(2,3,8,12)
#' pc2 <- pcor(base = th, p=4, parametrization='itp', itheta = ith)
#' pc2
#'
#' Sparse(solve(pc2$base), zeros.rm = TRUE)
#'
pcor.numeric <- function(base, p,
                         parametrization = 'cpc',
                         itheta, d0,
                         decomposition,
                         ...) {
  if(missing(p)) {
    if(missing(itheta)) {
      p <- (1 + sqrt(1+8*length(base)))/2
      stopifnot(floor(p)==ceiling(p))
    } else {
      stop("please provid 'p'!")
    }
  }
  stopifnot(p>1)
  if(missing(itheta)) {
    stopifnot(length(base)==(p*(p-1)/2))
    itheta <- which(lower.tri(diag(p)))
  }
  stopifnot(length(base)==length(itheta))
  theta <- base
  parametrization <- match.arg(
    arg = tolower(parametrization),
    choices = c("cpc", "sap", "itp")
  )
  if(missing(decomposition)) {
    decomposition <- "svd"
  }
  decomposition <- match.arg(
    arg = tolower(decomposition),
    choices = c("svd", "eigen", "chol")
  )
  if(missing(d0)) {
    d0 <- p:1
  }
  L <- theta2L(
    theta = theta,
    p = p,
    parametrization = parametrization,
    itheta = itheta,
    d0 = d0
  )
  if(attr(L, "parametrization") == 'itp') {
    C0 <- cov2cor(chol2inv(t(L)))
  } else {
    C0 <- tcrossprod(L)
  }
  return(
    pcor(
      base = C0,
      p = p,
      parametrization = parametrization,
      itheta = itheta,
      d0 = d0,
      decomposition = decomposition,
      ...)
    )
}
#' @rdname pcor-class
#' @export
pcor.matrix <- function(base, p,
                        parametrization = 'cpc',
                        itheta, d0,
                        decomposition,
                        ...) {
  parametrization <- match.arg(
    arg = tolower(parametrization),
    choices = c("cpc", "sap", "itp")
  )
  stopifnot(all.equal(base, t(base)))
  p <- as.integer(nrow(base))
  if(missing(itheta)) {
    itheta <- which(lower.tri(diag(p)))
  }
  if(missing(d0)) {
    d0 <- p:1
  }
  if(missing(decomposition)) {
    decomposition <- "svd"
  }
  decomposition <- match.arg(
    arg = tolower(decomposition),
    choices = c("svd", "eigen", "chol")
  )
  if(parametrization == 'itp') {
    Q <- chol2inv(chol(base))
    ilQ <- intersect(
      which(lower.tri(matrix(1, p, p))),
      which(!is.zero(Q))
    )
    if(missing(itheta)) {
      itheta <- ilQ
    } else {
      stopifnot(all(itheta %in% ilQ))
    }
    L <- t(chol(Q))
    if(missing(d0)) {
      d0 <- diag(L)
    } else {
      stopifnot((length(d0)==p) && (all(d0>0)))
      for(i in 1:p) {
        L[i, ] <- L[i, ]/L[i, i]
        L[i, ] <- L[i, ]*d0[i]
      }
    }
    theta <- L[itheta]
  } else {
    il <- which(lower.tri(base))
    l <- t(chol(base))[il]
    theta <- optim(
      rep(0, length(il)),
      function(x)
        mean((theta2L(x, p, parametrization)[il]-l)^2),
      method = 'BFGS')$par
  }
  out <- list(
    base = base,
    theta = theta,
    p = p,
    parametrization = parametrization,
    itheta = itheta,
    d0 = d0)
  out$H <- Hcorrel(
    theta = theta,
    p = p,
    parametrization = parametrization,
    itheta = itheta,
    d0 = d0,
    C0 = base,
    decomposition = decomposition,
    ...)
  class(out) <- "pcor"
  return(out)
}
#' @describeIn pcor
#' Print method for 'pcor'
#' @export
print.pcor <- function(x, ...) {
  cat("Parameters (", toupper(x$parametrization),
      " parametrization):\n", sep = "")
  print(x$theta)
  cat("Base correlation matrix:\n")
  print(x$base)
}
