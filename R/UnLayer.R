UnLayer <- function(snet, color_pos = "green3", color_neg = "red3") {
  if("static.sign" %in% class(snet)) {
    multi <- snet

    pos <- multi[["gal"]][[".subnetcache"]][[".LayerID"]][["+"]]
    neg <- multi[["gal"]][[".subnetcache"]][[".LayerID"]][["-"]]

    mat <- as.sociomatrix(pos) + as.sociomatrix(neg)

    net <- as.network(mat,matrix.type = "adjacency", directed = pos$gal$directed, loops = pos$gal$loops)

     #add edge attributes
     for (e in list.edge.attributes(multi)) {
       net%e%e <- multi%e%e
     }

    #add edge color
    net%e%'col' <- ifelse(net%e%'sign'== 1, color_pos, color_neg)

     #add vertex attributes
     for (v in list.vertex.attributes(pos)) {
       net%v%v <- pos%v%v
     }

    #save multilayer
    net$mlt <- multi

     comb <- net
  } else if ("dynamic.sign" %in% class(snet))  {
    multi <- c(snet[["gal"]][[".subnetcache"]][[".NetworkID"]][[1]][["gal"]][[".PrevNets"]],
                    snet[["gal"]][[".subnetcache"]][[".NetworkID"]])
    comb <- lapply(multi, function(x) {
      pos <- x[["gal"]][[".subnetcache"]][[".LayerID"]][["+"]]
      neg <- x[["gal"]][[".subnetcache"]][[".LayerID"]][["-"]]

      mat <- as.sociomatrix(pos) + as.sociomatrix(neg)

      net <- as.network(mat,matrix.type = "adjacency", directed = pos$gal$directed, loops = pos$gal$loops)

      #add edge attributes
      for (e in list.edge.attributes(snet)) {
        net%e%e <- x%e%e
      }

      #add edge color
      net%e%'col' <- ifelse(net%e%'sign'== 1, color_pos, color_neg)

      #add vertex attributes
      for (v in list.vertex.attributes(pos)) {
        net%v%v <- pos%v%v
      }

      #save multilayer
      net$mlt <- x

      x <- net
    })


  } else{
    stop("The input should be a network created by the signNetwork() function")
  }
  return(comb)
}
