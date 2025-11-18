#' ergm.sign: A Package for Exponential Random Graph Models for Signed Networks
#'
#' The ergm.sign package implements tools to simulate and estimate Signed Exponential Random Graph Models
#' and Temporal Signed Exponential Random Graph Models.
#'
#' @author
#' Marc Schalberger
#'
#' @import network
#' @importFrom ergm.multi Layer uncombine_network subnetwork_templates Networks
#' @importFrom ergm ergm ergmMPLE summary_formula ergm.getnetwork
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom tergm NetSeries control.tergm
#' @importFrom Rdpack reprompt
#' @importFrom intergraph asIgraph
#' @importFrom graphlayouts layout_with_stress
#' @importFrom vegan procrustes
#' @importFrom igraph induced_subgraph E E<-
#' @importFrom purrr map
#' @importFrom graphics boxplot lines
#' @importFrom methods is
#' @importFrom utils getFromNamespace packageVersion
#'
#' @name ergm.sign
NULL

