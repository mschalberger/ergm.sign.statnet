UnLayer <- function(snet) {
  if("static.sign" %in% class(snet)) {
    multi <- snet

    pos <- multi[["gal"]][[".subnetcache"]][[".LayerID"]][["+"]]
    neg <- multi[["gal"]][[".subnetcache"]][[".LayerID"]][["-"]]



    mat <- as.sociomatrix(pos, attrname = "sign") + as.sociomatrix(neg, attrname = "sign")

     net <- as.network(abs(mat),matrix.type = "adjacency", directed = pos$gal$directed, loops = pos$gal$loops)
     net <- set.edge.value(net, "sign", mat)

     comb <- net
    #comb <- sum(pos,neg, na.rm = T)
  } else if ("dynamic.sign" %in% class(snet))  {
    multi <- c(snet[["gal"]][[".subnetcache"]][[".NetworkID"]][[1]][["gal"]][[".PrevNets"]],
                    snet[["gal"]][[".subnetcache"]][[".NetworkID"]])
    comb <- lapply(multi_nets, function(x) {
      pos <- x[["gal"]][[".subnetcache"]][[".LayerID"]][["+"]]
      neg <- x[["gal"]][[".subnetcache"]][[".LayerID"]][["-"]]

      mat <- as.sociomatrix(pos, attrname = "sign") + as.sociomatrix(neg, attrname = "sign")

      net <- as.network(abs(mat),matrix.type = "adjacency", directed = pos$gal$directed, loops = pos$gal$loops)
      net <- set.edge.value(net, "sign", mat)

      x <- net
      #comb <- sum(pos, neg, na.rm = T)
    })
  } else{
    stop("The input should be a network created by the signNetwork() function")
  }
  return(list(
    multi = multi,
    single = comb))
}
