#' Germany data on four disease with geometry
#' @usage data(Germany4diseases)
#' @format sf with an 544 areas and 9 variables
#' @description
#' This is an `sf` object containing the observed
#' and expected number of cases on oral, oesophagus,
#' larynx and lung.
#' @keywords disease mapping, cancer, Germany
#' @examples
#' data(Germany4diseases)
#' ggplot2::ggplot(Germany4diseases) +
#'   ggplot2::geom_sf(aes(fill = lung_obs/lung_exp))
