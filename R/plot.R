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
#' @param net A signed network object of class \code{static.sign}.
#' @param col_pos Color for positive edges. Default is 'green3'.
#' @param col_neg Color for negative edges. Default is 'red3'.
#' @param neg.lty Line type for negative edges. Default is "dotted". Other options are "solid" and "dashed".
#' @param inv_weights Logical. If TRUE, edge weights are inverted (1/weights) so positive edges pull nodes closer together. Default is TRUE.
#' @param ... Additional arguments passed to the plot function.
#'
#' @return A plot of the signed network.
#'
#' @references
#' Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization.
#' *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#'
#' @seealso \link{UnLayer}, \link[graphlayouts]{layout_with_stress}
#'
#' @examples
#' data("tribes")
#' plot.static.sign(tribes, col_pos = "green", col_neg = "red")
#'
#' @export
plot.static.sign <- function(net, col_pos = "green3", col_neg = "red3",
                             neg.lty = 2, inv_weights = TRUE, ...) {
  sgl <- UnLayer(net, color_pos = col_pos, color_neg = col_neg, neg.lty = neg.lty)

  tmp_graph <- intergraph::asIgraph(sgl)
  weights <- if (inv_weights) 1 / sgl%e%"weights" else sgl%e%"weights"
  layout <- graphlayouts::layout_with_stress(tmp_graph, weights = weights)

  p <- network::plot.network(sgl,
            edge.col = "col",
            edge.lty = "type",
            coord = layout,
            ...)
  return(p)
}

#' Visualization for Dynamic Signed Networks
#'
#' \code{plot.dynamic.sign()} visualizes a dynamic signed network over multiple timepoints.
#'
#' @inheritParams plot.static.sign
#' @inheritSection plot.static.sign Layout
#' @param net A signed network object of class \code{dynamic.sign}.
#' @param time A vector of integers indicating which timepoints should be visualized. Defaults to all.
#' @param titles A character vector of names for the timepoints.
#' @param fix.pos Logical. If TRUE, the layout is fixed across timepoints based on the first timepoint. Default is TRUE.
#'
#' @return A list of plots, one for each selected timepoint.
#'
#' @export
plot.dynamic.sign <- function(net, col_pos = "#008000", col_neg = "#E3000F",
                              neg.lty = 2, inv_weights = TRUE,
                              time = NULL, titles = NULL, fix.pos = TRUE, ...) {
  sgl <- UnLayer(net, color_pos = col_pos, color_neg = col_neg, neg.lty = neg.lty)

  if (is.null(time)) {
    time <- seq_len(length(sgl))
  }

  p <- vector("list", length(time))

  if (fix.pos) {
    # Reference layout: t = 1
    ref_sgl <- sgl[[1]]
    ref_graph <- intergraph::asIgraph(ref_sgl)
    weights_ref <- if (inv_weights) 1 / ref_sgl%e%"weights" else ref_sgl%e%"weights"
    E(ref_graph)$weights <- weights_ref

    layout_ref <- graphlayouts::layout_with_stress(ref_graph, weights = E(ref_graph)$weights)

    for (t in time) {
      tmp_graph <- intergraph::asIgraph(sgl[[t]])
      tmp_weights <- if (inv_weights) 1 / sgl[[t]]%e%"weights" else sgl[[t]]%e%"weights"
      E(tmp_graph)$weights <- tmp_weights

      # Full stress layout for current graph
      layout_t <- graphlayouts::layout_with_stress(tmp_graph, weights = E(tmp_graph)$weights)

      # Rotate/translate the entire layout to match reference
      Yrot <- vegan::procrustes(layout_ref, layout_t, scale = F)$Yrot
      layout_t <- Yrot

      # Plot using rotated layout
      p[[t]] <- plot(sgl[[t]],
                     edge.col = "col",
                     edge.lty = "type",
                     coord = layout_t,
                     main = if (!is.null(titles)) titles[t] else NULL,
                     ...)
    }


  } else {
    p <- vector("list", length(time))
    for (t in time) {
      ref_sgl <- sgl[[t]]
      ref_graph <- intergraph::asIgraph(ref_sgl)
      weights <- if (inv_weights) 1 / ref_sgl%e%"weights" else ref_sgl%e%"weights"
      layout <- graphlayouts::layout_with_stress(ref_graph, weights = weights)
      p[[t]] <- plot(sgl[[t]],
                     edge.col = "col",
                     edge.lty = "type",
                     coord = layout,
                     main = titles[t],
                     ...)
    }
  }
  return(p)
}


