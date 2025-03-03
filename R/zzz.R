#  File R/zzz.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' @import ergm
#' @import network
#' @import statnet.common
#' @import stats
#' @importFrom Rdpack reprompt
## #' @import rlang

.onAttach <- function(libname, pkgname){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm.multi", c("statnet"), FALSE)
  options(ergm.ABI.action = "disable")  # Ignore ABI mismatch warnings
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}

.onLoad <- function(libname, pkgname){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")
  options(ergm.eval.loglik=TRUE)

  .RegisterProposals()
  .RegisterKeywords()

  InitErgmTerm..layer.net <- getFromNamespace("InitErgmTerm..layer.net", ns="ergm.multi")
  body(InitErgmTerm..layer.net)[[7]] <- quote(list(name="_layer_net",
                                                   coef.names = c(),
                                                   iinputs = c(unlist(.block_vertexmap(nw, ".LayerID", FALSE)),
                                                               if(is.directed(nw)) sapply(nwl, function(nw) (nw %v% ".undirected")[1]), ll),
                                                   dependence = dependence))
  assignInNamespace("InitErgmTerm..layer.net", InitErgmTerm..layer.net, ns="ergm.multi")
}

## BEGIN boilerplate: should be kept in sync with statnet.common.
# TODO: Figure out some automatic way to keep this in sync with statnet.common.
#' @name snctrl
#'
#' @title Statnet Control
#'
#' @description A utility to facilitate argument completion of control lists, reexported from `statnet.common`.
#'
#' @section Currently recognised control parameters:
#' This list is updated as packages are loaded and unloaded.
#'
#' \Sexpr[results=rd,stage=render]{statnet.common::snctrl_names()}
#'
#' @seealso [statnet.common::snctrl()]
#' @docType import
NULL
#' @export
snctrl <- statnet.common::snctrl
## END boilerplate: should be kept in sync with statnet.common.

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", "|bd|bdmax|.dyads&fixL", 0, "random", "randomtoggleFixL")
  ergm_proposal_table("c", "Bernoulli", "|bd|bdmax|.dyads&fixL&sparse", 1, "TNT", "TNTFixL")
}


.RegisterKeywords <- function() {
  ergm_keyword(name="layer-aware", short="layer", description="operates on multilayer network constructs", popular=TRUE, package="ergm.multi")
}

#' @useDynLib ergm.sign
.onUnload <- function(libpath){
  library.dynam.unload("ergm.sign",libpath)
}
