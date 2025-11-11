#' @keywords internal
"_PACKAGE"

#' ergm.sign: A Package for Exponential Random Graph Models for Signed Networks
#'
#' The ergm.sign package implements tools to simulate and estimate Signed Exponential Random Graph Models
#' and Temporal Signed Exponential Random Graph Models.
#'
#' @author
#' Marc Schalberger
#'
#' @import network
#' @import ergm
#' @import ergm.multi
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom tergm NetSeries control.tergm
#' @importFrom Rdpack reprompt
#' @importFrom intergraph asIgraph
#' @importFrom graphlayouts layout_with_stress
#' @importFrom vegan procrustes
#' @importFrom igraph induced_subgraph E
#' @importFrom purrr map
#' @name ergm.sign
NULL

