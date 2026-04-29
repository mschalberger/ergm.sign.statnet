#' @export
simulate.sign <- function(object, nsim = 1, seed = NULL, ...) {
  net <- object$network

  raw_sims <- NextMethod()

  if(nsim == 1) {
    sim <- raw_sims
    nwl <- uncombine_network(sim)
    #nwl <- c(nwl[[1L]] %n% ".PrevNets", nwl)
    sim %n% "NetList" <- nwl
  } else {
    sim <- lapply(raw_sims, function(x) {
      nwl <- uncombine_network(x)
      #nwl <- c(nwl[[1L]] %n% ".PrevNets", nwl)
      x %n% "NetList" <- nwl
      return(x)
    })
  }

}
