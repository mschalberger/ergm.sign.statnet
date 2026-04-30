#' @export
simulate.sign <- function(object, nsim = 1, seed = NULL, ...) {
  net <- object$network
  is_dynamic <- "dynamic.sign" %in% class(net)

  raw_sims <- NextMethod()
if(is_dynamic) {
  if(nsim == 1) {
    sim <- raw_sims
    nwl <- uncombine_network(sim)
    nwl <- lapply(nwl, function(x) {class(x) <- c("static.sign",class(x)); return(x)})
    sim %n% "NetList" <- nwl
  } else {
    sim <- lapply(raw_sims, function(x) {
      nwl <- uncombine_network(x)
      nwl <- lapply(nwl, function(y) {class(y) <- c("static.sign",class(y)); return(y)})
      x %n% "NetList" <- nwl
      return(x)
    })
  }
  }else {
  sim <- raw_sims
}

}
