#' Network Attributes for Signed Networks
#'
#' Print descriptive statistics of a signed network.
#'
#' @param net A signed network object of class \code{static.sign} or \code{dynamic.sign}.
#' @param time A list of integers indicating what timepoints should be summarised.
#'
#' @return Matrix with network attributes.
#'
#' @seealso \link{signNetwork}
#'
#' @examples
#' data("tribes")
#' summary(tribes)
#'
#' @export
summary.dynamic.sign <- function(net, time = NULL, names = NULL) {
  nws <- net$gal$NetList
  if (is.null(time)) {
    time <- 1:length(nws)
  }

  mat <- data.frame()
  for (i in time) {
    nw <- UnLayer(nws[[i]])
    MultiNet <- nws[[i]]
    n <- network.size(nw)
    a <- data.frame(Directed = nw$gal[["directed"]],
                    Loops = nw$gal[["loops"]],
                    Nodes = network.size(nw),
                    Edges = summary_formula(MultiNet ~ edges),
                    "Edges+" = summary_formula(MultiNet ~ L(~ edges, ~ `+`)),
                    "Edges-" = summary_formula(MultiNet ~ L(~ edges, ~ `-`)),
                    Triads = summary_formula(MultiNet ~ L(~ triangle, ~ `+` | `-`)),
                    `+++` = summary_formula(MultiNet ~ L(~ triangle, ~ `+`)),
                    `---` = summary_formula(MultiNet ~ L(~ triangle, ~ `-`)),
                    `++-` = ifelse(nw$gal[["directed"]],
                                   sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OTP") +
                                                         espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ITP")+
                                                         espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OSP") +
                                                         espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ISP")) *c(1:n)),
                          sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * c(1:n))),
                    `+--` = ifelse(nw$gal[["directed"]],
                                   sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OTP") +
                                                         espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ITP")+
                                                         espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OSP") +
                                                         espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ISP")) *c(1:n)),
                                   sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * c(1:n))),
                    Density = round(network.density(nw), 2), check.names = F)
    rownames(a) <- ifelse(is.null(names), paste("Time", i), names[i])
    mat <- rbind(mat, a)
  }
  return(mat)
}


#' @export
summary.static.sign <- function(net) {
  net <- UnLayer(net)
  MultiNet <- net$gal$mlt

  a <- as.data.frame(c(as.character(net$gal[["directed"]]),
                       as.character(net$gal[["loops"]]),
                       n <- network.size(net),
                       summary_formula(MultiNet~edges +
                                         L(~edges, ~`+`) +
                                         L(~edges, ~`-`) +
                                         L(~triangle, ~ `+`|`-`) +
                                         L(~triangle,~ `+`) +
                                         L(~triangle, ~ `-`)),
                       ifelse(net$gal[["directed"]],
                                sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OTP") +
                                                      espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ITP")+
                                                      espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OSP") +
                                                      espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ISP")) *c(1:n)),
                                sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * c(1:n))),
                       ifelse(net$gal[["directed"]],
                              sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OTP") +
                                                    espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ITP")+
                                                    espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OSP") +
                                                    espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ISP")) *c(1:n)),
                              sum(summary_formula(MultiNet ~ espL(d = c(1:n), L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * c(1:n))),
                       round(network.density(net),2)))
  rownames(a) <- c("Directed", "Loops","Nodes","Edges","Edges+", "Edges-", "Triads", "   +++", "   ---", "   ++-","   +--", "Density")
  colnames(a) <- ""
  cat("Network Attributes:")
  return(a)
}

