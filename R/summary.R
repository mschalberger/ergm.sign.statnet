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
  sgl      <- UnLayer(object)
  A        <- as.sociomatrix(sgl, attrname = "sign")
  directed <- isTRUE(sgl %n% "directed")

  s   <- .summary_from_adj(A, directed)
  out <- .format_summary_row(s)
  rownames(out) <- ""
  cat("Network Attributes:\n")
  return(out)
}

#' @rdname summary.static.sign
#' @param object A signed network object of class \code{dynamic.sign}.
#' @param time Integer vector of timepoints to summarize. Defaults to all.
#' @export
summary.dynamic.sign <- function(object, time = NULL, ...) {
  nws <- object %n% "NetList"
  if (is.null(time)) time <- seq_along(nws)

  mat <- do.call(rbind, lapply(time, function(i) {
    sgl      <- UnLayer(nws[[i]])
    A        <- as.sociomatrix(sgl, attrname = "sign")
    directed <- isTRUE(sgl %n% "directed")
    s        <- .summary_from_adj(A, directed)
    .format_summary_row(s)
  }))

  rownames(mat) <- NULL
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

# =============================================================================
# Matrix-based summary helpers
# =============================================================================

#' @keywords internal
.summary_from_adj <- function(A, directed) {
  n  <- nrow(A)

  # signed layers
  Ap <- (A ==  1L | A == 2) * 1L
  An <- (A == -1L | A == 2) * 1L
  Ab <- (abs(A) > 0) * 1L

  if (directed) {

    # --- all 2-path types (OTP, ITP, OSP, ISP)
    .all_path_types <- function(M) {
      Mt  <- t(M)
      OTP <- M  %*% M
      ITP <- Mt %*% Mt
      OSP <- M  %*% Mt
      ISP <- Mt %*% M
      OTP + ITP + OSP + ISP
    }

    paths_all <- .all_path_types(Ab)
    tri_all   <- sum(paths_all * Ab) / 3

    PP <- .all_path_types(Ap)
    NN <- .all_path_types(An)

    tri_ppp <- sum(PP * Ap) / 3 # + + +
    tri_nnn <- sum(NN * An) / 3 # - - -

    tri_ppm <- sum(PP * An)  # + + closed by -
    tri_pmm <- sum(NN * Ap)  # - - closed by +

    ep <- sum(Ap)
    en <- sum(An)
    e  <- ep + en

  } else {
    Ap2 <- Ap %*% Ap
    An2 <- An %*% An
    Ab2 <- Ab %*% Ab

    tri_ppp <- sum(diag(Ap2 %*% Ap)) / 6
    tri_nnn <- sum(diag(An2 %*% An)) / 6

    tri_ppm <- sum(Ap2 * An) / 2
    tri_pmm <- sum(An2 * Ap) / 2

    tri_all <- sum(diag(Ab2 %*% Ab)) / 6

    ep <- sum(Ap) / 2
    en <- sum(An) / 2
    e  <- ep + en
  }

  # density
  possible <- if (directed) n * (n - 1L) else n * (n - 1L) / 2
  density  <- round(e / possible, 3)

  list(
    directed = directed,
    loops    = FALSE,
    n        = n,
    edges    = e,
    edges_p  = ep,
    edges_n  = en,
    tri_all  = tri_all,
    tri_ppp  = tri_ppp,
    tri_nnn  = tri_nnn,
    tri_ppm  = tri_ppm,
    tri_pmm  = tri_pmm,
    density  = density
  )
}

#' @keywords internal
.format_summary_row <- function(s) {
  data.frame(
    Directed  = s$directed,
    Loops     = s$loops,
    Nodes     = s$n,
    Edges     = s$edges,
    `Edges+`  = s$edges_p,
    `Edges-`  = s$edges_n,
    Triads    = s$tri_all,
    `+++`     = s$tri_ppp,
    `---`     = s$tri_nnn,
    `++-`     = s$tri_ppm,
    `+--`     = s$tri_pmm,
    Density   = s$density,
    check.names = FALSE
  )
}

