#' Construct a Signed Network Representation for ERGM
#'
#' @description
#' `Signed()` is a convenience wrapper around [Layer()] that constructs a
#' two-layer multilayer network representation of a signed network, with one
#' layer for positive ties (`pos`) and one for negative ties (`neg`). It
#' accepts several input formats and handles the layer construction and
#' constraint specification automatically.
#'
#' @param ... Input networks or edge attributes. Three call signatures are
#'   supported:
#'   \describe{
#'     \item{`Signed(pos_nw, neg_nw)`}{Two separate [`network`] objects, the
#'       first containing positive ties and the second negative ties.}
#'     \item{`Signed(nw, sign_attr)`}{A single [`network`] object and the name
#'       of an edge attribute (character string) whose values encode the sign of
#'       each tie. Values must be in \{-1, 0, 1\}, where `1` indicates a
#'       positive tie, `-1` a negative tie, and `0` (or `NA`) an unsigned or
#'       absent tie.}
#'     \item{`Signed(nw, pos_attr, neg_attr)`}{A single [`network`] object and
#'       the names of two separate binary edge attributes encoding positive and
#'       negative ties respectively.}
#'   }
#' @param dual.sign Logical. If `FALSE` (the default), a
#'   constraint is added to prevent any dyad from simultaneously having both a
#'   positive and a negative tie. Set to `TRUE` to allow dual-signed ties, e.g.
#'   in multiplex contexts where the same dyad may appear in both layers.
#' @param .symmetric Passed to [Layer()]. A logical vector indicating which
#'   layers should be treated as undirected (symmetrized) even if the
#'   underlying network is directed.
#' @param .bipartite Passed to [Layer()]. An integer vector indicating the
#'   bipartite block size for each layer.
#' @param .active Passed to [Layer()]. A list of vertex attribute
#'   specifications, one per layer, indicating which vertices are active in
#'   each layer.
#'
#' @return A combined multilayer [`network`] object as returned by [Layer()],
#'   with two layers named `"pos"` and `"neg"`.
#'
#' @seealso [Layer()] for the underlying multilayer construction;
#'   [`ergm`][ergm::ergm] for model fitting on the resulting network.
#'
#' @examples
#' # Two separate network objects
#' fit <- ergm.sign(Signed(pos_nw, neg_nw) ~ Pos(~edges) + Neg(~edges))
#'
#' # Single network with a {-1, 0, 1}-valued edge attribute
#' fit <- ergm.sign(Signed(nw, "sign") ~ Pos(~edges) + Neg(~edges))
#'
#' # Single network with two separate binary edge attributes
#' fit <- ergm.sign(Signed(nw, "is_pos", "is_neg") ~ Pos(~edges) + Neg(~edges))
#'
#' # Allow dual-signed ties (pos and neg simultaneously)
#' fit <- ergm.sign(Signed(nw, "sign", dual.sign = TRUE) ~ Pos(~edges) + Neg(~edges))
#' @export
Signed <- function(..., dual.sign = FALSE, .symmetric = NULL, .bipartite = NULL, .active = NULL) {
  args <- list(...)

  if (length(args) == 2 && all(sapply(args, is, "network"))) {
    pos_nw <- args[[1]]
    neg_nw <- args[[2]]
    nwl <- list(`+` = pos_nw, `-` = neg_nw)
    MultiNet <- Layer(nwl, .symmetric = .symmetric, .bipartite = .bipartite, .active = .active)
  }

  if (length(args) >= 1 && is(args[[1]], "network")) {
    nw <- args[[1]]

    if (length(args) == 3 && is.character(args[[2]]) && is.character(args[[3]])) {
      pos_attr <- args[[2]]
      neg_attr <- args[[3]]
      MultiNet <- Layer(nw, c(`+` = pos_attr, `-` = neg_attr),
                   .symmetric = .symmetric, .bipartite = .bipartite, .active = .active)
    }

    # Sub-case 2b: Single edge attribute with binary values (1 / -1 or 1 / 0)
    if (length(args) == 2 && is.character(args[[2]])) {
      sign_attr <- args[[2]]
      sign_vals <- get.edge.attribute(nw, sign_attr)

      if (is.null(sign_vals))
        stop("Edge attribute ", sQuote(sign_attr), " not found in network.")

      unique_vals <- sort(unique(sign_vals[!is.na(sign_vals)]))

      if (!all(unique_vals %in% c(-1, 0, 1)))
        stop(
          "Edge attribute ", sQuote(sign_attr),
          " must contain only values in {-1, 0, 1}. ",
          "For two separate attributes, pass them as two character arguments."
        )

      # Derive pos/neg binary edge attributes
      nw %e% "pos" <- as.integer(sign_vals == 1)
      nw %e% "neg" <- as.integer(sign_vals == -1)

      MultiNet <- Layer(nw, c(`+` = "pos", `-` = "neg"),
                        .symmetric = .symmetric, .bipartite = .bipartite, .active = .active)

    }
  }

  if (exists("MultiNet")) {
    if (!dual.sign) MultiNet %ergmlhs% "constraints" <- update(MultiNet %ergmlhs% "constraints", ~. + ChangeStats(~L(~edges, ~`+` & `-`)))
    MultiNet %v% "sign" <- MultiNet %v% ".LayerName"
    class(MultiNet) <- c("static.sign", "network", class(MultiNet))
    MultiNet%n%"dual.sign" <- dual.sign
    return(MultiNet)
  }

  stop(
    "Unrecognized format for Signed(). Valid calls:\n",
    "  Signed(pos_nw, neg_nw)                     # two network objects\n",
    "  Signed(nw, 'sign_attr')                     # one attr with values in {-1, 0, 1}\n",
    "  Signed(nw, 'pos_attr', 'neg_attr')          # two separate binary edge attributes"
  )
}
