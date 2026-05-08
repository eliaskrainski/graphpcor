#' @importFrom methods new getMethod
#' @importFrom grDevices rgb
#' @importFrom utils tail
#' @importFrom stats vcov cov2cor runif rexp as.formula drop1
#' @importFrom Matrix Matrix t forceSymmetric colSums Diagonal
#' @importFrom numDeriv hessian
#' @importFrom INLAtools is.zero Sparse upperPadding packageCheck
#' @importFrom INLAtools cgeneric cgeneric_Q cgeneric_graph
#' @importFrom INLAtools cgeneric_initial cgeneric_prior
#' @importFrom igraph graph_from_adjacency_matrix
#' @useDynLib graphpcor, .registration = TRUE
NULL
