#' Visualization for Signed Networks
#'
#' Functions to visualize signed networks in static or dynamic form.
#'
#' @section Layout:
#' Uses a force-directed graph layout based on stress majorization, implemented in the \code{graphlayouts} package
#' via \code{layout_with_stress()}. Similar to Kamada-Kawai, but generally faster and with better results.
#'
#' @section Static signed networks:
#' \code{plot.static.sign()} visualizes a single (static) signed network.
#'
#' @param x A signed network object of class \code{static.sign}.
#' @param col_pos Color for positive edges. Default is 'green3'.
#' @param col_neg Color for negative edges. Default is 'red3'.
#' @param neg.lty Line type for negative edges. Default is "solid". Other options are "dotted" and "dashed".
#' @param inv_weights Logical. If TRUE, edge weights are inverted (1/weights) so positive edges pull nodes closer together. Default is TRUE.
#' @param coord Optional matrix of coordinates for node positions. If NULL, layout is computed using stress majorization.
#' @param ... Additional arguments passed to the plot function.
#'
#' @return A plot of the signed network.
#'
#' @references \insertRef{gansner2004graph}{ergm.sign}
#'
#' @seealso \link{UnLayer}, \link[graphlayouts]{layout_with_stress}
#'
#' @examples
#' data("tribes")
#' plot(tribes, col_pos = "green", col_neg = "red")
#'
#' @export
plot.static.sign <- function(x, col_pos = "#008000", col_neg = "#E3000F",
                             neg.lty = 1, inv_weights = TRUE, coord = NULL, ...) {
  net <- x
  sgl <- UnLayer(net, color_pos = col_pos, color_neg = col_neg, neg.lty = neg.lty)

  tmp_graph <- intergraph::asIgraph(sgl)
  weights <- sgl%e%"weights"
  if (inv_weights) weights <- 1 / pmax(weights, .Machine$double.eps)

  # Only compute layout if no coord is provided
  if (is.null(coord)) {
    coord <- graphlayouts::layout_with_stress(tmp_graph, weights = weights)
  }

  network::plot.network(
    sgl,
    edge.col = sgl%e%"col",
    edge.lty = sgl%e%"type",
    coord = coord,
    ...
  )
}


#' Visualization for Dynamic Signed Networks
#'
#' \code{plot.dynamic.sign()} visualizes a dynamic signed network over multiple timepoints.
#'
#' @inheritParams plot.static.sign
#' @inheritSection plot.static.sign Layout
#' @param x A signed network object of class \code{dynamic.sign}.
#' @param time A vector of integers indicating which timepoints should be visualized. Defaults to all.
#' @param titles A character vector of names for the timepoints.
#' @param fix.pos Logical. If TRUE, the layout is fixed across timepoints based on the first timepoint. Default is TRUE.
#' @param ... Additional arguments passed to the plot function.
#'
#' @return A list of plots, one for each selected timepoint.
#'
#' @export
plot.dynamic.sign <- function(x, col_pos = "#008000", col_neg = "#E3000F",
                              neg.lty = 1, inv_weights = TRUE,
                              time = NULL, titles = NULL, fix.pos = TRUE, ...) {
  net <- x
  sgl <- UnLayer(net, color_pos = col_pos, color_neg = col_neg, neg.lty = neg.lty)

  if (is.null(time)) time <- seq_along(sgl)
  p <- vector("list", length(time))

  if (fix.pos) {
    # Reference layout from first timepoint
    ref_sgl <- sgl[[1]]
    ref_graph <- intergraph::asIgraph(ref_sgl)
    weights_ref <- ref_sgl%e%"weights"
    if (inv_weights) weights_ref <- 1 / pmax(weights_ref, .Machine$double.eps)
    E(ref_graph)$weights <- weights_ref

    layout_ref <- graphlayouts::layout_with_stress(ref_graph, weights = E(ref_graph)$weights)

    for (t in time) {
      tmp_graph <- intergraph::asIgraph(sgl[[t]])
      tmp_weights <- sgl[[t]]%e%"weights"
      if (inv_weights) tmp_weights <- 1 / pmax(tmp_weights, .Machine$double.eps)
      E(tmp_graph)$weights <- tmp_weights

      layout_t <- graphlayouts::layout_with_stress(tmp_graph, weights = E(tmp_graph)$weights)

      # Align layout to reference
      layout_t <- vegan::procrustes(layout_ref, layout_t, scale = FALSE)$Yrot

      p[[t]] <- network::plot.network(sgl[[t]],
                                      edge.col = sgl[[t]]%e%"col",
                                      edge.lty = sgl[[t]]%e%"type",
                                      coord = layout_t,
                                      main = if (!is.null(titles)) titles[t] else NULL,
                                      ...)
    }
  } else {
    for (t in time) {
      tmp_sgl <- sgl[[t]]
      tmp_graph <- intergraph::asIgraph(tmp_sgl)
      tmp_weights <- tmp_sgl%e%"weights"
      if (inv_weights) tmp_weights <- 1 / pmax(tmp_weights, .Machine$double.eps)
      E(tmp_graph)$weights <- tmp_weights

      layout <- graphlayouts::layout_with_stress(tmp_graph, weights = E(tmp_graph)$weights)

      p[[t]] <- network::plot.network(tmp_sgl,
                                      edge.col = tmp_sgl%e%"col",
                                      edge.lty = tmp_sgl%e%"type",
                                      coord = layout,
                                      main = if (!is.null(titles)) titles[t] else NULL,
                                      ...)
    }
  }

  p
}


