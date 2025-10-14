#' A signed multi- or dynamic-network representation.
#'
#' Combines multiple signed networks into a single object for ERGM modeling.
#' Can represent either multilayer or dynamic network structures.
#'
#' @param ... Several signed networks or a list of signed networks.
#' @param dynamic Logical. If \code{TRUE}, treat the input as a dynamic network; otherwise, as a multilayer network.
#' @return A combined network object with appropriate constraints for model fitting.
#' @export


signNetworks <- function(..., dynamic = FALSE) {
  args <- list(...)
  if(all(sapply(args, is, "static.sign"))){
    nwl <- as.list(args)
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "static.sign"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multinetwork specification. See help for information.")

  # cycle through time and add lagged attributes
  # for (i in seq_along(nwl)) {
  #   nw <- nwl[[i]]
  #
  #   # vertex attributes
  #   vattrs <- list.vertex.attributes(nw)
  #   if (i > 1) {
  #     prev <- nwl[[i-1]]
  #     for (nm in vattrs) {
  #       nw %v% paste0("lag_", nm) <- prev %v% nm
  #     }
  #   } else {
  #     for (nm in vattrs) {
  #       nw %v% paste0("lag_", nm) <- NA
  #     }
  #   }
  #
  #   # edge attributes
  #   eattrs <- list.edge.attributes(nw)
  #   if (i > 1) {
  #     prev <- nwl[[i-1]]
  #     for (nm in eattrs) {
  #       nw %e% paste0("lag_", nm) <- prev %e% nm
  #     }
  #   } else {
  #     for (nm in eattrs) {
  #       nw %e% paste0("lag_", nm) <- NA
  #     }
  #   }
  #
  #   nwl[[i]] <- nw
  # }
  #
  # nwl[[1]]%n%"lag_neg" <-  matrix(NA, nrow = network.size(nwl[[1]]), ncol = network.size(nwl[[1]]))  # first time step has no lagged neg
  # for (i in 2:length(nwl)) {
  #   prev <- nwl[[i-1]]
  #   mat <- as.sociomatrix(prev, attrname = "neg")  # numeric 0/1 matrix
  #   nwl[[i]]%n%"lag_neg"<- mat                  # network attribute for this time step
  # }

  # Remove constraints
  nwl <- lapply(nwl, function(nw) {
    nw%ergmlhs%"constraints" <- NULL
    return(nw)
  })

  # Combine networks
  comb <- if (dynamic) NetSeries(nwl) else Networks(nwl)
  group_var <- if (dynamic) ".TimeID" else ".NetworkID"

  #comb[["gal"]][[".subnetcache"]] <- NULL

  # Create grouping ID for blockdiag
  comb %v% ".NetworkID_new" <- as.numeric(factor(paste(comb %v% group_var, comb%v%".LayerID", sep = "-")))

  # Add constraint
  new_constraints <- if (dynamic) {
    as.formula("~ fixL(~`+` & `-`) + blockdiag('.NetworkID_new') + discord('.PrevNet')")
  } else {
    as.formula("~ fixL(~`+` & `-`) + blockdiag('.NetworkID_new')")
  }
  comb %ergmlhs% "constraints" <- new_constraints
  return(comb)
}
