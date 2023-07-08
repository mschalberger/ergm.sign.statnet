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
summary.dynamic.sign <- function(net, time = c(1:length(net)), names = NULL) {
  class(net) <- "network"
  mat <- data.frame()
  for (i in time) {
    nw <- net[[i]]
    n <- network.size(nw)
    MultiNet <- Layer(nw, c(`+` = "pos", `-` = "neg"))
    a <- data.frame(Directed = nw$gal[["directed"]],
                    Loops = nw$gal[["loops"]],
                    Nodes = network.size(nw),
                    Edges = summary_formula(MultiNet ~ edges),
                    `+ edges` = summary_formula(MultiNet ~ L(~ edges, ~ `+`)),
                    `- edges` = summary_formula(MultiNet ~ L(~ edges, ~ `-`)),
                    Triads = summary_formula(MultiNet ~ L(~ triangle, ~ `+` | `-`)),
                    `+++` = summary_formula(MultiNet ~ L(~ triangle, ~ `+`)),
                    `---` = summary_formula(MultiNet ~ L(~ triangle, ~ `-`)),
                    `++-` = sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * c(1:n)),
                    `+--` = sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * c(1:n)),
                    Density = round(network.density(nw), 2), check.names = F)
    rownames(a) <- ifelse(is.null(names), paste("Time", i), names[i])
    mat <- rbind(mat, a)
  }
  #rownames(mat) <- ifelse(is.null(names), time, names)
  return(mat)
}


#' @rdname summary
#' @export
summary.static.sign <- function(net) {
  class(net) <- "network"
  MultiNet <- Layer(net, c(`+` = "pos",`-`= "neg"))
  a <- as.data.frame(c(as.character(net$gal[["directed"]]),
                       as.character(net$gal[["loops"]]),
                       n <- network.size(net),
                       summary_formula(MultiNet~edges +
                                         L(~edges, ~`+`) +
                                         L(~edges, ~`-`) +
                                         L(~triangle, ~ `+`|`-`) +
                                         L(~triangle,~ `+`) +
                                         L(~triangle, ~ `-`)),
                       sum(summary_formula(MultiNet ~ espL(c(1:(n)), L.base= ~`-`, Ls.path= c(~`+`,~`+`))) * c(1:(n))),
                       sum(summary_formula(MultiNet ~ espL(d = c(1:(n)), L.base= ~`+`, Ls.path= c(~`-`,~`-`))) * c(1:(n))),
                       round(network.density(net),2)))
  rownames(a) <- c("Directed", "Loops","Nodes","Edges","   + edges", "   - edges", "Triads", "   +++", "   ---", "   ++-","   +--", "Density")
  colnames(a) <- ""
  cat("Network Attributes:")
  return(a)
}

