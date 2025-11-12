#' @templateVar name randomtoggleFixL
#' @aliases InitErgmProposal.randomtoggleFixL
#' @title Propose a randomly selected dyad to toggle, respecting the layer constraint
#' @description Propose a randomly selected dyad to toggle
#' @template ergmProposal-general
NULL
InitErgmProposal.randomtoggleFixL <- function(arguments, nw){
  set_layer_namemap <- utils::getFromNamespace(".set_layer_namemap", "ergm.multi")
  mk_layer_net_auxform <- utils::getFromNamespace(".mk_.layer.net_auxform", "ergm.multi")

  Ls <- set_layer_namemap(arguments$constraints$fixL$Ls, nw)
  if(is(Ls, "formula")) Ls <- list(Ls)
  auxiliaries <- mk_layer_net_auxform(Ls)

  list(name = "randomtoggleFixL", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw), auxiliaries = auxiliaries)
}

#' @templateVar name TNTFixL
#' @aliases InitErgmProposal.TNTFixL
#' @title Default MH algorithm respecting the layer constraint
#' @description Stratifies the population of dyads
#'   edge status: those having ties and those having no ties (hence T/NT).
#'   This is useful for improving performance in sparse networks,
#'   because it gives at least 50\% chance of proposing a toggle of an existing edge.
#' @template ergmProposal-general
NULL
InitErgmProposal.TNTFixL <- function(nw, arguments, ...){
  set_layer_namemap <- utils::getFromNamespace(".set_layer_namemap", "ergm.multi")
  mk_layer_net_auxform <- utils::getFromNamespace(".mk_.layer.net_auxform", "ergm.multi")

  Ls <- set_layer_namemap(arguments$constraints$fixL$Ls, nw)
  if(is(Ls, "formula")) Ls <- list(Ls)
  auxiliaries <- mk_layer_net_auxform(Ls)


  list(name = "TNTFixL", dyadgen = ergm_dyadgen_select(arguments, nw), bd = ergm_bd_init(arguments, nw), auxiliaries = auxiliaries)
}
