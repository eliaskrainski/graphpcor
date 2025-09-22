#' Set a tree whose nodes represent the two kind of variables:
#' children and parent.
setClass("treepcor")

#' Set a graph whose nodes and edges represent variables and
#' conditional distributions, respectively.
setClass("graphpcor")

#' Information for a base model for correlation matrices
setClass(
  "basecor",
  slots = c("base", "theta", "p", "parametrization",
            "itheta", "H"),
  validity = function(object) {
    p3 <- c("cpc", "CPC", "sap", "SAP", "itp", "ITP")
    (object$p>1) &&
      (object$p == nrow(object$base)) &&
      all.equal(object$base == t(object$base)) &&
      all(diag(object$base==1)) &&
      all(diag(chol(object$base))>0) &&
      object$parametrization %in% p3
  }
)
