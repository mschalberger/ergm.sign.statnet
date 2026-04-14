#' Network Attributes for Signed Networks
#'
#' Print descriptive statistics of a signed network.
#'
#' @section Static signed networks:
#' \code{summary.static.sign()} summarizes a single (static) signed network.
#'
#' @param object A signed network object of class \code{static.sign}.
#' @param ... Additional arguments.
#'
#' @return A data frame or matrix with network attributes.
#'
#' @seealso \link{network.sign}, \link{UnLayer}
#'
#' @examples
#' data("tribes")
#' summary(tribes)
#'
#' @export
summary.static.sign <- function(object, ...) {
  net <- object
  net_sgl <- UnLayer(net)
  multi <- net_sgl%n%"mlt"
  n <- network.size(net_sgl)

  a <- data.frame(
    Directed = net_sgl%n%"directed",
    Loops = net_sgl%n%"loops",
    Nodes = n,
    Edges = summary_formula(multi ~ edges),
    `Edges+` = summary_formula(multi ~ L(~ edges, ~ `+`)),
    `Edges-` = summary_formula(multi ~ L(~ edges, ~ `-`)),
    Triads = summary_formula(multi ~ L(~ triangle, ~ `+` | `-`)),
    `+++` = summary_formula(multi ~ L(~ triangle, ~ `+`)),
    `---` = summary_formula(multi ~ L(~ triangle, ~ `-`)),
    `++-` = if (net_sgl%n%"directed") {
      sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OTP") +
                            espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ITP") +
                            espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OSP") +
                            espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ISP")) * (1:n))
    } else {
      sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * (1:n))
    },
    `+--` = if (net_sgl%n%"directed") {
      sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OTP") +
                            espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ITP") +
                            espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OSP") +
                            espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ISP")) * (1:n))
    } else {
      sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * (1:n))
    },
    Density = round(network.density(net_sgl), 2),
    check.names = FALSE
  )
  rownames(a) <- ""
  cat("Network Attributes:\n")
  return(a)
}

#' @rdname summary.static.sign
#' @param object A signed network object of class \code{dynamic.sign}.
#' @param time Integer vector of timepoints to summarize. Defaults to all.
#' @export
summary.dynamic.sign <- function(object, time = NULL, ...) {
  net <- object
  nws <- net%n%"NetList"
  if (is.null(time)) time <- seq_along(nws)

  mat <- do.call(rbind, lapply(time, function(i) {
    nw_sgl <- UnLayer(nws[[i]])
    multi <- nws[[i]]
    n <- network.size(nw_sgl)

    row <- data.frame(
      Directed = nw_sgl%n%"directed",
      Loops = nw_sgl%n%"loops",
      Nodes = n,
      Edges = summary_formula(multi ~ edges),
      `Edges+` = summary_formula(multi ~ L(~ edges, ~ `+`)),
      `Edges-` = summary_formula(multi ~ L(~ edges, ~ `-`)),
      Triads = summary_formula(multi ~ L(~ triangle, ~ `+` | `-`)),
      `+++` = summary_formula(multi ~ L(~ triangle, ~ `+`)),
      `---` = summary_formula(multi ~ L(~ triangle, ~ `-`)),
      `++-` = if (nw_sgl%n%"directed") {
        sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OTP") +
                              espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ITP") +
                              espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "OSP") +
                              espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`), type = "ISP")) * (1:n))
      } else {
        sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `-`, Ls.path = c(~ `+`, ~ `+`))) * (1:n))
      },
      `+--` = if (nw_sgl%n%"directed") {
        sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OTP") +
                              espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ITP") +
                              espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "OSP") +
                              espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`), type = "ISP")) * (1:n))
      } else {
        sum(summary_formula(multi ~ espL(d = 1:n, L.base = ~ `+`, Ls.path = c(~ `-`, ~ `-`))) * (1:n))
      },
      Density = round(network.density(nw_sgl), 2),
      check.names = FALSE
    )
    rownames(row) <- NULL
    row
  }))

  return(mat)
}

#' Summary formula method for dynamic signed networks
#'
#' Calculates statistics for dynamic.sign objects at specified timepoints.
#'
#' @param object A formula with a dynamic.sign network as LHS.
#' @param at Numeric vector of timepoints. Defaults to all if missing.
#' @param basis Optional dynamic.sign network. If NULL, uses LHS network.
#' @param ... Additional arguments passed to summary_formula for network objects.
#' @return Matrix of statistics for each timepoint.
#' @importFrom ergm summary_formula
#' @importFrom utils getS3method
#' @export
summary_formula.dynamic.sign <- function(object, at,  ..., basis = NULL) {
  basis <- if (is.null(basis)) ergm.getnetwork(object) else basis
  if (missing(at) || !is.numeric(at)) at <- seq_along(basis%n%"NetList")

  res_list <- lapply(at, function(t) {
    nw <- basis%n%"NetList"
    nw <- nw[[t]]
    getS3method("summary_formula", "network")(object, basis = nw, ...)
  })

  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  return(res)
}

