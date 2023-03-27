#' Network Attributes for Signed Networks
#'
#' Print descriptive statistics of a signed network.
#'
#' @param net A signed network object of class \code{static.sign} or \code{dynamic.sign}.
#' @param time A list of integers indicating what timepoints should be summarised.
#'
#' @return Matrix with network attributes.
#'
#' @seealso \link{signnet}
#'
#' @export

summary <- function(net, time = c(1:length(net))) UseMethod("summary")

#' @rdname summary
#' @export
summary.dynamic.sign <- function(net, time = c(1:length(net))) {
  class(net) <- "network"
    mat <- c()
    for (i in time) {
      nw <- net[[i]]
      n = network.size(nw)
      MultiNet <- Layer(nw, c(`+` = "pos",`-`= "neg"))
      a <- matrix(c(as.character(nw$gal[["directed"]]),
                    nw$gal[["loops"]],
                    network.size(nw),
                    summary_formula(MultiNet~edges +
                                      L(~edges, ~`+`) +
                                      L(~edges, ~`-`) +
                                      L(~triangle, ~ `+`|`-`) +
                                      L(~triangle,~ `+`) +
                                      L(~triangle, ~ `-`)),
                    sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base= ~`-`, Ls.path= c(~`+`,~`+`))) * c(1:n)),
                    sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base= ~`+`, Ls.path= c(~`-`,~`-`))) * c(1:n)),
                    round(network.density(nw),2)), ncol =1)
      rownames(a) <- c("Directed", "Loops", "Multiple","Nodes","Edges","   + edges", "   - edges", "Triads", "   +++", "   ---", "   ++-","   +--", "Density")
      colnames(a) <- paste("Time", i)
      mat <- cbind(mat,a)
    }
    print(noquote(mat))
}

#' @rdname summary
#' @export
summary.static.sign <- function(net) {
  class(net) <- "network"
  MultiNet <- Layer(net, c(`+` = "pos",`-`= "neg"))
  a <- matrix(c(as.character(net$gal[["directed"]]),
                net$gal[["loops"]],
                n <- network.size(net),
                summary_formula(MultiNet~edges +
                                  L(~edges, ~`+`) +
                                  L(~edges, ~`-`) +
                                  L(~triangle, ~ `+`|`-`) +
                                  L(~triangle,~ `+`) +
                                  L(~triangle, ~ `-`)),
                sum(summary_formula(MultiNet ~ espL(c(1:(n)), L.base= ~`-`, Ls.path= c(~`+`,~`+`))) * c(1:(n))),
                sum(summary_formula(MultiNet ~ espL(d = c(1:(n)), L.base= ~`+`, Ls.path= c(~`-`,~`-`))) * c(1:(n))),
                round(network.density(net),2)), ncol =1)
  rownames(a) <- c("Directed", "Loops", "Multiple","Nodes","Edges","   + edges", "   - edges", "Triads", "   +++", "   ---", "   ++-","   +--", "Density")
  colnames(a) <- ""
  cat("Network Attributes:")
  return(noquote(a))
}
