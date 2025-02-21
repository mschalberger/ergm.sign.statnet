#'  A signed multinetwork network representation.
#'
#'  A function for specifying the LHS of a signed multi-network (a.k.a. multilevel) ERGM. Typically used in conjunction with the \code{N()} term operator.
#'
#'  @param ... Several signed networks or a list of signed networks.
#'
#'  @return A network object comprising the provided networks, with multinetwork metadata.
#'
#'  @export

signNetworks <- function(...) {
  args <- list(...)
  if(all(sapply(args, is, "static.sign"))){
    nwl <- as.list(args)
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "static.sign"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multinetwork specification. See help for information.")

  # Remove layer constraint
  nwl <- lapply(nwl, function(nw) nw%ergmlh%"constraints" <- formula(~fixL(~`+` & `-`)))

  # Combine networks
  comb <- Networks(nwl)

  # Add new blockdiag attribute
  comb%v%".NetworkID_new" <-  as.numeric(factor(paste0(comb%v%".LayerID", comb%v%".NetworkID")))

  # Add constraint
  comb%ergmlh%"constraints" <- formula(~fixL(~`+` & `-`) + blockdiag(".NetworkID_new"))

  return(comb)
}
