#' Multilayer network to single layer network.
#'
#' Turn a multilayer network object into a single layer network object.
#'
#' @param net A signed network object of class \code{static.sign} or \code{dynamic.sign}.
#' @param color_pos Color for positive edges. Default is '#008000'.
#' @param color_neg Color for negative edges. Default is '#E3000F'.
#' @param color_both Color for edges that are both positive and negative. Default is '#333333'.
#' @param neg.lty Line type for negative edges. Default is 2.
#' @param both.lty Line type for edges that are both positive and negative. Default is 2.
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
UnLayer <- function(net, color_pos = "#008000", color_neg = "#E3000F", color_both= "#333333", neg.lty = 2, both.lty = 2) {

  if ("static.sign" %in% class(net)) {

    pos <- get.inducedSubgraph(net, v = which(net %v% ".LayerName" == "+"))
    neg <- get.inducedSubgraph(net, v = which(net %v% ".LayerName" == "-"))

    mat_pos <- as.sociomatrix(pos)
    mat_neg <- as.sociomatrix(neg)

    adj <- ((mat_pos + mat_neg) > 0) * 1

    comb <- as.network(adj, matrix.type = "adjacency",
                       directed = pos$gal$directed, loops = pos$gal$loops)

    # Copy edge attributes
    for (e in list.edge.attributes(net)) {
      comb%e%e <- net%e%e
    }

    # Sign attribute: 1 = pos only, -1 = neg only, 2 = both
    sign_attr <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
    sign_attr[mat_pos > 0 & mat_neg == 0] <- 1
    sign_attr[mat_pos == 0 & mat_neg > 0] <- -1
    sign_attr[mat_pos > 0 & mat_neg > 0] <- 2

    comb%e%'sign' <- sign_attr
    comb%e%'col' <- ifelse(sign_attr == 1, color_pos,
                              ifelse(sign_attr == -1, color_neg, color_both))
    comb%e%'weights' <- ifelse(sign_attr == 1, 2,
                                  ifelse(sign_attr == -1, 1, 1.5))
    comb%e%'type' <- ifelse(sign_attr == 1, 1,
                               ifelse(sign_attr == -1, neg.lty, both.lty))

    attributes <- setdiff(list.vertex.attributes(pos), c(".LayerName", ".LayerID"))

    # Copy vertex attributes
    for (v in attributes) {
      comb%v%v <- pos%v%v
    }

    # Store original multilayer network
    comb$gal$mlt <- net

  } else if (any(c("dynamic.sign", "multi.sign") %in% class(net))) {

    multi <- net$gal$NetList
    comb <- lapply(multi, function(x) {

      pos <- get.inducedSubgraph(x, v = which(x %v% ".LayerName" == "+"))
      neg <- get.inducedSubgraph(x, v = which(x %v% ".LayerName" == "-"))

      mat_pos <- as.sociomatrix(pos)
      mat_neg <- as.sociomatrix(neg)

      adj <- ((mat_pos + mat_neg) > 0) * 1

      net_sgl <- as.network(adj, matrix.type = "adjacency",
                            directed = pos$gal$directed, loops = pos$gal$loops)

      # Copy edge attributes
      for (e in list.edge.attributes(x)) {
        net_sgl%e%e <- x%e%e
      }

      sign_attr <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
      sign_attr[mat_pos > 0 & mat_neg == 0] <- 1
      sign_attr[mat_pos == 0 & mat_neg > 0] <- -1
      sign_attr[mat_pos > 0 & mat_neg > 0] <- 2

      net_sgl%e%'sign' <- sign_attr
      net_sgl%e%'col' <- ifelse(sign_attr == 1, color_pos,
                                ifelse(sign_attr == -1, color_neg, color_both))
      net_sgl%e%'weights' <- ifelse(sign_attr == 1, 2,
                                    ifelse(sign_attr == -1, 1, 1.5))
      net_sgl%e%'type' <- ifelse(sign_attr == 1, 1,
                                 ifelse(sign_attr == -1, neg.lty, both.lty))

      attributes <- setdiff(list.vertex.attributes(pos), c(".LayerName", ".LayerID"))

      # Copy vertex attributes
      for (v in attributes) {
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

