InitErgmTerm.Pos <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  formula <- if(is(a$formula, "formula")) list(a$formula) else a$formula

  a$formula <- if(length(formula)==1) formula[[1]] else formula
  cl <- call("L", a$formula, Ls = ~`+`)
  trm <- call.ErgmTerm(cl, nw, ...)

  trm
}

InitErgmTerm.Neg <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  formula <- if(is(a$formula, "formula")) list(a$formula) else a$formula

  a$formula <- if(length(formula)==1) formula[[1]] else formula
  cl <- call("L", a$formula, Ls = ~`-`)
  trm <- call.ErgmTerm(cl, nw, ...)

  trm
}

#' @templateVar name dsf
#' @title Dyadwise shared friends
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of dyads in the network with exactly `d[i]` shared friends. For a directed network, multiple shared friends definitions are possible.
#'
#' @usage
#' # binary: dsf(d, type="OTP", in_order=FALSE)
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.dsf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","in_order"),
                      vartypes = c("numeric","character","logical"),
                      defaultvalues = list(NULL, "OTP", FALSE),
                      required = c(TRUE, FALSE,FALSE))
  cl <- call("ddspL", d = a$d, type = a$type, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names<- paste("dsf", a$d, sep = ".")
  trm
}

#' @templateVar name dse
#' @title Dyadwise shared enemies
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of dyads in the network with exactly `d[i]` shared enemies. For a directed network, multiple shared enemies definitions are possible.
#'
#' @usage
#' # binary: dse(d, type="OTP", in_order=FALSE)
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.dse <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","in_order"),
                      vartypes = c("numeric","character","logical"),
                      defaultvalues = list(NULL, "OTP", FALSE),
                      required = c(TRUE, FALSE,FALSE))
  cl <- call("ddspL", d = a$d, type = a$type, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- paste("dse", a$d, sep = ".")
  trm
}

#' @templateVar name gwdsf
#' @title Geometrically weighted dyadwise shared friends distribution
#' @description This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared friends distribution with decay parameter. Note that the GWDSF statistic is equal to the sum of GWNSF plus GWESF. For a directed network, multiple shared friend definitions are possible.
#'
#' @usage
#' # binary: gwdsf(decay, fixed=FALSE, cutoff=30, type="OTP", in_order=FALSE)
#' @templateVar multiplicand shared friend or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying DSF
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.gwdsf <- function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP",NULL,FALSE),
                      required = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))
  cl <- call("dgwdspL", decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type, alpha = a$alpha, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  #trm$coef.names <- paste("gwdsf", a$decay, sep = ".")
  trm
}

#' @templateVar name gwdse
#' @title Geometrically weighted dyadwise shared enemies distribution
#' @description This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared enemies distribution with decay parameter. Note that the GWDSE statistic is equal to the sum of GWNSE plus GWESE. For a directed network, multiple shared friend definitions are possible.
#'
#' @usage
#' # binary: gwdse(decay, fixed=FALSE, cutoff=30, type="OTP", in_order=FALSE)
#' @templateVar multiplicand shared enemy or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying DSE
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.gwdse <- function(nw, arglist,cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP",NULL,FALSE),
                      required = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))
  cl <- call("dgwdspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type, alpha = a$alpha, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  #tmp$coef.names <- paste("gwdse", decay, sep = ".")
  trm
}

#' @templateVar name esf
#' @title Edgewise shared friends
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of edges in the network with exactly `d[i]` shared friends. For a directed network, multiple shared friend definitions are possible.
#'
#' @usage
#' # binary: esf(d, type="OTP", L.base=NULL, in_order=FALSE)
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.esf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "formula","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  cl <- call("despL", d = a$d, type = a$type, L.base = a$base, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name ese
#' @title Edgewise shared enemies
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of edges in the network with exactly `d[i]` shared enemies. For a directed network, multiple shared enemy definitions are possible.
#'
#' @usage
#' # binary: ese(d, type="OTP", L.base=NULL, in_order=FALSE)
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.esf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "formula","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  cl <- call("despL",d = a$d, type = a$type, L.base = a$base, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name gwesf
#' @title Geometrically weighted edgewise shared friend distribution
#' @description This term adds a statistic equal to the geometrically weighted edgewise (not dyadwise) shared friend distribution with decay parameter. For a directed network, multiple shared friend definitions are possible.
#'
#' @usage
#' # binary: gwesf(decay, fixed=FALSE, cutoff=30, type="OTP", base=NULL, in_order=FALSE)
#' @templateVar multiplicand shared friend or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying ESF
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.gwesf <- function(nw, arglist, cache.sp=TRUE, gw.cutoff=30,...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "formula","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  cl <- call("dnspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = a$base, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name gwese
#' @title Geometrically weighted edgewise shared enemy distribution
#' @description This term adds a statistic equal to the geometrically weighted edgewise (not dyadwise) shared enemy distribution with decay parameter. For a directed network, multiple shared enemy definitions are possible.
#'
#' @usage
#' # binary: gwese(decay, fixed=FALSE, cutoff=30, type="OTP", base=NULL, in_order=FALSE)
#' @templateVar multiplicand shared enemy or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying ESE
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.gwese <- function(nw, arglist, cache.sp=TRUE, gw.cutoff=30,...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "formula","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  cl <- call("dnspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = a$base, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name nsf
#' @title Non-edgewise shared friends
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of non-edges in the network with exactly `d[i]` shared friends. For a directed network, multiple shared friend definitions are possible.
#'
#' @usage
#' # binary: nsf(d, type="OTP", base=NULL, in_order=FALSE)
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.nsf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "formula","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  cl <- call("dnspL",d = a$d, type = a$type, L.base = a$base, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name nse
#' @title Non-edgewise shared enemies
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of non-edges in the network with exactly `d[i]` shared enemies. For a directed network, multiple shared enemy definitions are possible.
#'
#' @usage
#' # binary: nse(d, type="OTP", base=NULL, in_order=FALSE)
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.nse <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "formula","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  cl <- call("dnspL",d = a$d, type = a$type, L.base = a$base, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name gwnsf
#' @title Geometrically weighted non-edgewise shared friend distribution
#' @description This term adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared friend distribution with decay parameter. For a directed network, multiple shared friend definitions are possible.
#'
#' @usage
#' # binary: gwnsf(decay, fixed=FALSE, cutoff=30, type="OTP", base=NULL, in_order=FALSE)
#' @templateVar multiplicand shared friend or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying NSF
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.gwnsf <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "formula","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  cl <- call("dgwnspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = a$base, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

#' @templateVar name gwnse
#' @title Geometrically weighted non-edgewise shared enemey distribution
#' @description This term adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared enemy distribution with decay parameter. For a directed network, multiple shared enemy definitions are possible.
#'
#' @usage
#' # binary: gwnse(decay, fixed=FALSE, cutoff=30, type="OTP", base=NULL, in_order=FALSE)
#' @templateVar multiplicand shared enemy or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying NSE
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#' @template ergmTerm-L-base
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept layer-aware
InitErgmTerm.gwnse <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "formula","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  cl <- call("dgwnspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = a$base, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm
}

