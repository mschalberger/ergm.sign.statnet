#  File R/InitErgmConstraint.R in package ergm.multi, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

##########################################################################################
# Each of the <InitErgmConstraint.X> functions accepts an existing constraint list, 'conlist',
# and to this adds an empty constraint list for term X; if any arguments are passed besides
# 'conlist", execution will halt.
#
# --PARAMETERS--
#   conlist: a list, presumably of constraints for other terms
#
# --RETURNED--
#   conlist: updated to include the initialized empty constraint list for term X
#
##########################################################################################

#### This constraint has been provisionally moved to 'ergm'. It may be
#### moved back here once 'ergm.multi' is on CRAN.
## #' @templateVar name blockdiag
## #' @title Block-diagonal structure constraint
## #' @description Force a block-diagonal structure (and its bipartite analogue) on
## #'   the network. Only dyads \eqn{(i,j)} for which
## #'   `attr(i)==attr(j)` can have edges.
## #'
## #'   Note that the current implementation requires that blocks be
## #'   contiguous for unipartite graphs, and for bipartite
## #'   graphs, they must be contiguous within a partition and must have
## #'   the same ordering in both partitions. (They do not, however,
## #'   require that all blocks be represented in both partitions, but
## #'   those that overlap must have the same order.)
## #'
## #'   If multiple block-diagonal constraints are given, or if
## #'   `attr` is a vector with multiple attribute names, blocks
## #'   will be constructed on all attributes matching.
## #'
## #' @usage
## #' # blockdiag(attr)
## #' @template ergmTerm-attr
## #'
## #' @template ergmConstraint-general
## #'
## #' @concept dyad-independent
## #' @concept directed
## #' @concept undirected
## #' @import rle
## InitErgmConstraint.blockdiag<-function(lhs.nw, attr=NULL, ...){
##   if(length(list(...)))
##     stop(paste("Block diagonal constraint takes one argument at this time."), call.=FALSE)
##   list(attr=attr,
##        free_dyads = {
##          n <- network.size(lhs.nw)
##          storage.mode(n) <- "integer"
##          a <- c(ergm_get_vattr(attr, lhs.nw)) # Strip attributes, which confuse rle().
##          if(NVL(lhs.nw%n%"bipartite",0)){
##            bip <- lhs.nw %n% "bipartite"
##            ea <- a[seq_len(bip)]
##            aa <- a[bip+seq_len(n-bip)]
##            if(length(rle(ea)$lengths)!=length(unique(rle(ea)$values)) || length(rle(aa)$lengths)!=length(unique(rle(aa)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks of the egos and the alters be contiguous. See ", sQuote("ergmConstraint?blockdiag"), " for more information.")

##            tmp <- .double.rle(ea, aa)
##            el <- tmp$lengths1
##            al <- tmp$lengths2

##            o <- rlebdm(c(rep(rle(FALSE), bip*n, scale="run"),
##                          do.call(c,rep(
##                                      mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
##                                             el, cumsum(el), SIMPLIFY=FALSE),
##                                      al)
##                                  )), n)
##            # Future-proofing: in case it's bipartite directed, add
##            # both thte blocks and their transposes. (If undirected,
##            # it'll get filtered out by the .attributes constraints.)
##            ot <- rlebdm(c(do.call(c,rep(
##                                       mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bip+bend-blen, blen, n-bip-bend), scale="run")},
##                                              al, cumsum(al), SIMPLIFY=FALSE),
##                                       el)
##                                   ),
##                           rep(rle(FALSE), (n-bip)*n, scale="run")), n)
##            compress(o | ot)
##          }else{
##            a <- rle(a)
##            rlebdm(compress(do.call(c,rep(
##                                        mapply(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
##                                               a$lengths, cumsum(a$lengths), SIMPLIFY=FALSE),
##                                        a$lengths)
##                                    )), n)
##          }
##        },
##        dependence = FALSE)
## }

#' @templateVar name fixL
#' @title Logical layer constraint
#' @description This layer-aware constraint limits the sample space to those networks for which the specified logical layers are unchanged
#'
#' @usage
#' # fixL(Ls)
#' @templateVar Ls.interp , specifying fixed relations
## #' @template ergmTerm-Ls-1
#'
#' @template ergmConstraint-general
#'
#' @concept layer-aware
#' @concept directed
#' @concept undirected
InitErgmConstraint.fixL<-function(nw, arglist,  ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("Ls"),
                      vartypes = c("formula,list"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  list(Ls = a$Ls, dependence=TRUE)
}

