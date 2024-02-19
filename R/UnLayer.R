UnLayer <- function(snet) {
  if("static.sign" %in% class(snet)) {
    pos <- snet[["gal"]][[".subnetcache"]][[".LayerID"]][["+"]]
    neg <- snet[["gal"]][[".subnetcache"]][[".LayerID"]][["-"]]
    comb <- sum(pos,neg, na.rm = T)
  } else if ("dynamic.sign" %in% class(snet))  {
    multi_nets <- c(snet[["gal"]][[".subnetcache"]][[".NetworkID"]][[1]][["gal"]][[".PrevNets"]],
                    snet[["gal"]][[".subnetcache"]][[".NetworkID"]])
    comb <- lapply(multi_nets, function(x) {
      pos <- x[["gal"]][[".subnetcache"]][[".LayerID"]][["+"]]
      neg <- x[["gal"]][[".subnetcache"]][[".LayerID"]][["-"]]
      comb <- sum(pos, neg, na.rm = T)
    })
  } else{
    stop("The input should be a network created by the signNetwork() function")
  }
  return(comb)
}
