#' Network Attributes for Signed Networks
#'
#' Print descriptive statistics of a signed network.
#'
#' @section Static signed networks:
#' \code{summary.static.sign()} summarizes a single (static) signed network.
#'
#' @param net A signed network object of class \code{static.sign}.
#'
#' @return A data frame or matrix with network attributes.
#'
#' @seealso \link{network.sign}
#'
#' @examples
#' data("tribes")
#' summary.static.sign(tribes)
#'
#' @export
summary.static.sign <- function(net) {
  net <- UnLayer(net)
  MultiNet <- net$gal$mlt
  n <- network.size(net)

  a <- as.data.frame(t(c(
    as.character(net$gal[["directed"]]),
    as.character(net$gal[["loops"]]),
    n,
    summary_formula(MultiNet ~ edges +
                      L(~ edges, ~ `+`) +
                      L(~ edges, ~ `-`) +
                      L(~ triangle, ~ `+` | `-`) +
                      L(~ triangle, ~ `+`) +
                      L(~ triangle, ~ `-`)),
    if (net$gal[["directed"]]) {
      sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OTP") +
                            espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ITP") +
                            espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OSP") +
                            espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ISP")) * (1:n))
    } else {
      sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * (1:n))
    },
    if (net$gal[["directed"]]) {
      sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OTP") +
                            espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ITP") +
                            espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OSP") +
                            espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ISP")) * (1:n))
    } else {
      sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * (1:n))
    },
    round(network.density(net), 2)
  )))

  colnames(a) <- c("Directed", "Loops", "Nodes", "Edges", "Edges+", "Edges-",
                   "Triads", "+++", "---", "++-", " +--", "Density")
  rownames(a) <- ""
  cat("Network Attributes:\n")
  return(a)
}

#' @rdname summary.static.sign
#' @inheritParams summary.static.sign
#' @param time A vector of integers indicating which timepoints should be summarized.
#' @param names A character vector of names for the timepoints. If \code{NULL}, timepoints are numbered.
#' @export
summary.dynamic.sign <- function(net, time = NULL, names = NULL) {
  nws <- net$gal$NetList
  if (is.null(time)) time <- seq_along(nws)

  mat <- data.frame()
  for (i in time) {
    nw <- UnLayer(nws[[i]])
    MultiNet <- nws[[i]]
    n <- network.size(nw)

    a <- data.frame(
      Directed = nw$gal[["directed"]],
      Loops = nw$gal[["loops"]],
      Nodes = n,
      Edges = summary_formula(MultiNet ~ edges),
      "Edges+" = summary_formula(MultiNet ~ L(~ edges, ~ `+`)),
      "Edges-" = summary_formula(MultiNet ~ L(~ edges, ~ `-`)),
      Triads = summary_formula(MultiNet ~ L(~ triangle, ~ `+` | `-`)),
      `+++` = summary_formula(MultiNet ~ L(~ triangle, ~ `+`)),
      `---` = summary_formula(MultiNet ~ L(~ triangle, ~ `-`)),
      `++-` = if (nw$gal[["directed"]]) {
        sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OTP") +
                              espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ITP") +
                              espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OSP") +
                              espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ISP")) * (1:n))
      } else {
        sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * (1:n))
      },
      `+--` = if (nw$gal[["directed"]]) {
        sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OTP") +
                              espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ITP") +
                              espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OSP") +
                              espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ISP")) * (1:n))
      } else {
        sum(summary_formula(MultiNet ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * (1:n))
      },
      Density = round(network.density(nw), 2),
      check.names = FALSE
    )

    rownames(a) <- if (is.null(names)) paste("Time", i) else names[i]
    mat <- rbind(mat, a)
  }
  return(mat)
}

#' Calculation of summary statistics for dynamic signed netowrks
#'
#' A method for [summary_formula()] to calculate the
#' specified statistics for an observed [`dynamic.sign`] at the
#' specified time point(s).
#'
#' @aliases summary.formula
#' @param object An [`formula`] object with a
#' [`dynamic.sign`] as its LHS.
#' @param at Numeric vector specifying the time point(s) at which to
#' calculate the statistics. If missing, statistics are calculated at all time points.
#' @param basis An optional [`dynamic.sign`] object. If `NULL`, the network from the formula's LHS is used.
#' @param ... Additional arguments passed to or used by methods.
#' @return A matrix with \code{length(at)} row, one for each time
#' point in \code{at}, and columns for each term of the formula,
#' containing the corresponding statistics measured on the network
#' @importFrom ergm summary_formula
#' @importFrom utils getS3method
#' @export
summary_formula.dynamic.sign <- function(object, at, ..., basis=NULL) {
  basis <- NVL(basis, ergm.getnetwork(object))
  if (missing(at) || !is.numeric(at)) {
    at <- seq_along(basis$gal$NetList)
  }
  do.call(rbind, lapply(at, function(t) {
    nw <- basis$gal$NetList[[t]]
    getS3method("summary_formula","network")(object, basis = nw, ...)
  }))
}
