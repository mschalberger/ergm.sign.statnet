#' Multilayer network to single layer network.
#'
#' Turn a multilayer network object into a single layer network object.
#'
#' @param net A signed network object of class \code{static.sign} or \code{dynamic.sign}.
#' @param color_pos Color for positive edges. Default is '#008000'.
#' @param color_neg Color for negative edges. Default is '#E3000F'.
#' @param neg.lty Line type for negative edges. Default is 2.
#'
#' @return Single layer network object or a list of network objects for dynamic.sign.
#'
#' @seealso \link{network.sign}
#'
#' @examples
#' data("tribes")
#' tribes_sgl <- UnLayer(tribes)
#'
#' @export
UnLayer <- function(net, color_pos = "#008000", color_neg = "#E3000F", neg.lty = 2) {

  if ("static.sign" %in% class(net)) {

    pos <- get.inducedSubgraph(net, v = which(net %v% ".LayerName" == "+"))
    neg <- get.inducedSubgraph(net, v = which(net %v% ".LayerName" == "-"))

    mat <- as.sociomatrix(pos) - as.sociomatrix(neg)

    comb <- as.network(abs(mat), matrix.type = "adjacency",
                       directed = pos$gal$directed, loops = pos$gal$loops)

    # Copy edge attributes
    for (e in list.edge.attributes(net)) {
      comb%e%e <- net%e%e
    }

    # Add sign and visual attributes
    comb%e%'sign' <- mat
    comb%e%'col' <- ifelse(comb%e%'sign' == 1, color_pos, color_neg)
    comb%e%'weights' <- ifelse(comb%e%'sign' == 1, 2, 1)
    comb%e%'type' <- ifelse(comb%e%'sign' == 1, 1, neg.lty)

    # Copy vertex attributes
    for (v in list.vertex.attributes(pos)) {
      comb%v%v <- pos%v%v
    }

    # Store original multilayer network
    comb$gal$mlt <- net

  } else if ("dynamic.sign" %in% class(net)) {

    multi <- net$gal$NetList
    comb <- lapply(multi, function(x) {

      pos <- get.inducedSubgraph(x, v = which(x %v% ".LayerName" == "+"))
      neg <- get.inducedSubgraph(x, v = which(x %v% ".LayerName" == "-"))

      mat <- as.sociomatrix(pos) - as.sociomatrix(neg)

      net_sgl <- as.network(abs(mat), matrix.type = "adjacency",
                            directed = pos$gal$directed, loops = pos$gal$loops)

      # Copy edge attributes
      for (e in list.edge.attributes(x)) {
        net_sgl%e%e <- x%e%e
      }

      net_sgl%e%'sign' <- mat
      net_sgl%e%'col' <- ifelse(net_sgl%e%'sign' == 1, color_pos, color_neg)
      net_sgl%e%'weights' <- ifelse(net_sgl%e%'sign' == 1, 2, 1)
      net_sgl%e%'type' <- ifelse(net_sgl%e%'sign' == 1, 1, neg.lty)

      # Copy vertex attributes
      for (v in list.vertex.attributes(pos)) {
        net_sgl%v%v <- pos%v%v
      }

      # Store original multilayer network
      net_sgl$gal$mlt <- x

      net_sgl
    })

  } else {
    stop("Input must be a network created by network.sign()")
  }

  return(comb)
}

