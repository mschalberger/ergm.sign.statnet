#' Visualization for Signed Networks
#'
#' This function visualizes a signed network.
#'
#' @param net A signed network object of class \code{static.sign} or \code{dynamic.sign}.
#' @param col_pos Color for positive edges. Default is 'green3'.
#' @param col_neg Color for negative edges. Default is 'red3'.
#' @param neg.lyt Line type for negative edges. Default is "dotted". Other options are "solid" and "dashed".
#' @param inv_weights Logical. If TRUE, edge weights are inverted (1/weights) to ensure that the endpoints of positive edges are closer togehter than negative. Default is TRUE.
#' @param titles A character vector of names for the timepoints.
#' @param time A list of integers indicating what timepoints should be visualized.
#' @param ... Additional arguments passed to the plot function.
#' @return A plot of the signed network.
#' @seealso \link{UnLayer}
#' @examples
#' data("tribes")
#' plot.static.sign(tribes, col_pos = "green", col_neg = "red")
#'
#' @export
plot.static.sign <- function(net, col_pos = "green3", col_neg = "red3", neg.lyt = "dotted", inv_weights = T, ...) {
  sgl <- UnLayer(net, color_pos = col_pos, color_neg = col_neg, neg.lyt = neg.lyt)

  #get layout
  tmp_graph <- intergraph::asIgraph(sgl)
  weights <- if (inv_weights) {
    1 / sgl%e%"weights"
  } else {
    sgl%e%"weights"
  }
  layout <- graphlayouts::layout_with_stress(tmp_graph, weights = weights)

  #plot
  p <- plot(sgl,
       edge.col = "col",
       edge.lty = "type",
       coord = layout,
       ...)

  return(p)
}

#' @export
plot.dynamic.sign <- function(net, col_pos = "green3", col_neg = "red3", neg.lyt = "dotted", inv_weights = T, time = NULL, titles = NULL, fix.pos = T, ...) {
  sgl <- UnLayer(net, color_pos = col_pos, color_neg = col_neg, neg.lyt = neg.lyt)

  if (is.null(time)) {
    time <- seq_len(length(sgl))
  }

  if (fix.pos) {
    ref_sgl <- sgl[[1]]
    ref_graph <- intergraph::asIgraph(ref_sgl)

    # Compute edge weights
    weights <- if (inv_weights) {
      1 / ref_sgl%e%"weights"
    } else {
      ref_sgl%e%"weights"
    }

    # Identify connected nodes and isolates
    deg <- igraph::degree(ref_graph)
    connected <- which(deg > 0)
    isolates <- which(deg == 0)

    # Layout for connected nodes
    layout_conn <- graphlayouts::layout_with_stress(
      igraph::induced_subgraph(ref_graph, connected),
      weights = weights
    )

    # Compute bounding box of connected layout
    center <- colMeans(layout_conn)
    radius <- max(sqrt(rowSums((layout_conn - matrix(center, nrow(layout_conn), 2, byrow=TRUE))^2)))

    # Place isolates evenly spaced on a larger circle
    n_iso <- length(isolates)
    if (n_iso > 0) {
      angles <- seq(0, 2*pi, length.out = n_iso+1)[- (n_iso+1)]
      layout_iso <- cbind(center[1] + (radius * 1.2) * cos(angles),
                          center[2] + (radius * 1.2) * sin(angles))
    }

    # Combine layouts
    layout <- matrix(0, nrow = vcount(ref_graph), ncol = 2)
    layout[connected, ] <- layout_conn
    if (n_iso > 0) layout[isolates, ] <- layout_iso

    # Plot for each time step
    p <- vector("list", length(time))
    for (t in time) {
      p[[t]] <- plot(sgl[[t]],
                     edge.col = "col",
                     edge.lty = "type",
                     coord = layout,
                     main = titles[t],
                     ...)
    }
  } else {
    p <- c()
    for (t in time) {
      ref_sgl <- sgl[[t]]
      ref_graph <- intergraph::asIgraph(ref_sgl)
      weights <- if (inv_weights) {
        1 / ref_sgl%e%"weights"
      } else {
        ref_sgl%e%"weights"
      }
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
