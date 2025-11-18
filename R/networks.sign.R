#' Combine Signed Networks into a Multi- or Dynamic-Network Object
#'
#' Creates a composite network object from multiple signed networks,
#' suitable for ERGM modeling. Can represent either a multilayer
#' or dynamic signed network structure.
#'
#' @param ... One or more signed networks (objects of class \code{"static.sign"}),
#'   or a list of such networks.
#' @param dynamic Logical. If \code{TRUE}, treat input as a dynamic network;
#'   otherwise as a multilayer network. Defaults to \code{FALSE}.
#' @param dual.sign Logical. If \code{TRUE}, disables the layer fixing constraint. Defaults to \code{FALSE}.
#'
#' @return A combined network object of class \code{"multi.sign"} or
#'   \code{"dynamic.sign"}, with the appropriate ERGM constraint formula.
#'
#' @examples
#' data("tribes")
#' multi_net <- networks.sign(tribes, tribes)
#' dyn_net <- networks.sign(list(tribes, tribes), dynamic = TRUE)
#'
#' @export
networks.sign <- function(..., dynamic = FALSE, dual.sign = FALSE) {
  args <- list(...)

  # Validate input
  if (all(sapply(args, inherits, "static.sign"))) {
    nwl <- as.list(args)
  } else if (is.list(args[[1]]) && all(sapply(args[[1]], inherits, "static.sign"))) {
    nwl <- args[[1]]
  } else {
    stop("Invalid input format for multi-network specification. Expected signed network objects or a list of them.")
  }

  # Remove constraints from each network
  nwl <- lapply(nwl, function(nw) {
    nw %ergmlhs% "constraints" <- NULL
    nw
  })

  # Combine networks
  comb <- if (dynamic) NetSeries(nwl) else Networks(nwl)
  group_var <- if (dynamic) ".TimeID" else ".NetworkID"

  # Create grouping variable for blockdiag
  comb %v% ".NetworkID_new" <- as.numeric(factor(paste(
    comb %v% group_var,
    comb %v% ".LayerID",
    sep = "-"
  )))
  comb %v% ".LayerID" <- comb %v% ".NetworkID_new"

  ids <- unique(comb %v% ".LayerID")
  pairs <- split(ids, ceiling(seq_along(ids) / 2))

  # Build constraint expressions
  pair_str <- paste0(
    sapply(pairs, function(p) paste0("(`", p[1], "` & `", p[2], "`)")),
    collapse = " + "
  )
  fixL_str <- paste0("~ fixL(~", pair_str, ") + ")

  # Define ERGM constraints
  new_constraints <- if (dynamic) {
    as.formula(paste0(
      ifelse(dual.sign, "~", fixL_str),
      " blockdiag('.NetworkID_new') + discord('.PrevNet')"
    ))
  } else {
    as.formula(paste0(
      ifelse(dual.sign, "", fixL_str),
      " blockdiag('.NetworkID_new')"
    ))
  }

  comb %ergmlhs% "constraints" <- new_constraints

  # Assign network type class
  if (dynamic) {
    comb$gal$NetList <- nwl
    class(comb) <- c("dynamic.sign", class(comb))
  } else {
    class(comb) <- c("multi.sign", class(comb))
  }

  comb
}


