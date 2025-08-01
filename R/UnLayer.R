#' Multilayer network to single layer network.
#'
#' Turn a multilayer network object into a single layer network object.
#'
#' @param net A signed network object of class \code{static.sign} or \code{dynamic.sign}.
#' @param color_pos Color for positive edges. Default is 'green3'.
#' @param color_neg Color for negative edges. Default is 'red3'.
#'
#' @return Single layer network object.
#'
#' @seealso \link{signNetwork}
#'
#' @examples
#' data("tribes")
#' tribes_sgl <- UnLayer(tribes)
#'
#' @export
UnLayer <- function(net, color_pos = "green3", color_neg = "red3") {
  if("static.sign" %in% class(net)) {
    multi <- net

    pos <- get.inducedSubgraph(x = multi, v = which(multi%v%".LayerName"=="+"))
    neg <- get.inducedSubgraph(x = multi, v = which(multi%v%".LayerName"=="-"))

    mat <- as.sociomatrix(pos) + as.sociomatrix(neg) *-1

    net <- as.network(abs(mat),matrix.type = "adjacency", directed = pos$gal$directed, loops = pos$gal$loops)

     #add edge attributes
     for (e in list.edge.attributes(multi)) {
       net%e%e <- multi%e%e
     }

    net <- set.edge.value(net,'sign',mat)

    #add edge color
    net%e%'col' <- ifelse(net%e%'sign'== 1, color_pos, color_neg)

     #add vertex attributes
     for (v in list.vertex.attributes(pos)) {
       net%v%v <- pos%v%v
     }

    #save multilayer
    net$gal$mlt <- multi

     comb <- net
  } else if ("dynamic.sign" %in% class(net))  {
    multi <- net$gal$NetList
    comb <- lapply(multi, function(x) {
      pos <- get.inducedSubgraph(x = x, v = which(x%v%".LayerName"=="+"))
      neg <- get.inducedSubgraph(x = x, v = which(x%v%".LayerName"=="-"))

      mat <- as.sociomatrix(pos) + as.sociomatrix(neg) *-1

      net <- as.network(abs(mat),matrix.type = "adjacency", directed = pos$gal$directed, loops = pos$gal$loops)

      #add edge attributes
      for (e in list.edge.attributes(net)) {
        net%e%e <- x%e%e
      }

      if (!("sign" %in% list.edge.attributes(net))) {
        net%e%'sign' <- mat
      }


      #add edge color
      net%e%'col' <- ifelse(net%e%'sign'== 1, color_pos, color_neg)

      #add vertex attributes
      for (v in list.vertex.attributes(pos)) {
        net%v%v <- pos%v%v
      }

      #save multilayer
      net$gal$mlt <- x

      x <- net
    })


  } else{
    stop("The input should be a network created by the signNetwork() function")
  }
  return(comb)
}
