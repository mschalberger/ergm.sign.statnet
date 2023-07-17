#' Plot Signed Network for Signed Networks
#'
#' The function \code{plot.static.sign} or \code{plot.dynamic.sign} produces one or multiple simple two-dimensional plot of a dynamic signed network. In addition to the arguments of the \code{\link{network::plot.network}} function, some arguments for the context of signed networks have been added.
#'
#' @param net An object of class \code{static.sign} or \code{dynamic.sign} created with \code{\link{signNetwork}}.
#' @param ... Additional arguments to plot \code{\link{network::plot.network}}
#' @param time A vector of integers indicating which timepoints should be plotted.
#' @param color_pos Edge color for positive edges (1). Default is set to green.
#' @param color_neg Edge color for negative edges (-1). Default is set to red.
#' @param vertex.legend Display vertex legend, if \code{vertex.col} is a vertex attribute.
#' @param vertex.col Choose a vertex attribute as a color for the vertices.
#' @param vertex.legend.pos Specifying legend’s position using one of the following options: ”bottomright”, ”bottom”, ”bottomleft”, ”left”, ”topleft”, ”top”, ”topright”, ”right” or ”center”.
#' @param main A list of strings indicating the plot titles. If not specified, the plots will be named ”Time 1” to ”Time T”. If no titles are desired, the input should be set to NULL.
#' @param vertex.legend.size Size of the legend, which takes a character expansion factor relative to the current \code{par("cex")}.
#'
#' @seealso \link{signnet}
#'
#' @return  One or multiple two-dimensional plots.
#'
#' @export

plot <- function(net,
                 time = c(1:length(net)),
                 color_pos = "green3",
                 color_neg = "red3",
                 vertex.col = 2,
                 vertex.legend = F,
                 vertex.legend.pos = "topleft",
                 main = paste("Time ", c(1:length(net))),
                 vertex.legend.size = 0.65,
                 ...) UseMethod("plot")

#' @rdname plot
#' @export
plot.dynamic.sign <- function(net,
                             time = c(1:length(net)),
                             color_pos = "green3",
                             color_neg = "red3",
                             vertex.col = 2,
                             vertex.legend = F,
                             vertex.legend.pos = "topleft",
                             main = paste("Time ", c(1:length(net))),
                             vertex.legend.size = 0.65,
                             ...) {
    for (i in time) {
      nw <- net[[i]]
      nw%e%'sign' <- ifelse(nw%e%'sign'== 1, color_pos, color_neg)
      plot.network(nw,
                   main = main[i],
                   edge.col = "sign",
                   label = network.vertex.names(nw),
                   vertex.col = vertex.col,
                   ... = ...)
      if (vertex.legend == TRUE) {
        unique <- table(get.vertex.attribute(nw, vertex.col))
        legend(vertex.legend.pos, fill = c(0:length(unique)+1),
               legend = names(unique), cex = vertex.legend.size )
      }
    }
}

#' @rdname plot
#' @export
plot.static.sign <- function(net,
                             color_pos = "green3",
                             color_neg = "red3",
                             vertex.legend = F,
                             vertex.col = 2,
                             vertex.legend.pos = "topleft",
                             legend.size = 0.65,
                             ...) {
  class(net) <- "network"
  net%e%'sign' <- ifelse(net%e%'sign'== 1, color_pos, color_neg)
  plot.network(net,
               edge.col = "sign",
               label = network.vertex.names(net),
               vertex.col = vertex.col,
               ... = ...)
  if (vertex.legend == TRUE) {
    unique <- table(get.vertex.attribute(net, vertex.col))
    legend(vertex.legend.pos, fill = c(0:length(unique)+1),
           legend = names(unique), cex = legend.size)
  }
}

