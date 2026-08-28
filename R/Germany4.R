#' Germany data on four disease with geometry
#' @name Germany4
#' @usage data(Germany4)
#' @format sf with an 544 areas and 9 variables
#' @description
#' This is an `sf` object containing the observed
#' and expected number of cases on oral, oesophagus,
#' larynx and lung.
#' @docType data
#' @keywords dataset
#' @keywords Germany
#' @keywords cancer
#' @keywords disease
#' @keywords mapping
#' @keywords map
#' @examples
#' data(Germany4)
#' names(Germany4)
#' \dontrun{
#' ## ggplot2::ggplot(Germany4) +
#'   ## ggplot2::geom_sf(aes(fill = lung_obs/lung_exp))
#' }
NULL
"Germany4"
#' @rdname Germany4
#' @aliases graphGermany
"graphGermany"
