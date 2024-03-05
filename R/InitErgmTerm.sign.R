#' @templateVar name Pos
#' @title Evaluation of positive edges
#' @description Evaluates the terms in `formula` of the positive edges and sums the results elementwise.
#'
#' @usage
#' # binary: Pos(formula)
#'
#' @template ergmTerm-formula
#'
#' @concept operator
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
  trm$coef.names <- gsub("^.*~", "Pos~", trm$coef.names)
  trm
}

#' @templateVar name Neg
#' @title Evaluation of negative edges
#' @description Evaluates the terms in `formula` of the negative edges and sums the results elementwise.
#'
#' @usage
#' # binary: Neg(formula)
#'
#' @template ergmTerm-formula
#'
#' @concept operator
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
  trm$coef.names <- gsub("^.*~", "Neg~", trm$coef.names)
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
InitErgmTerm.dsf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","in_order"),
                      vartypes = c("numeric","character","logical"),
                      defaultvalues = list(NULL, "OTP", FALSE),
                      required = c(TRUE, FALSE,FALSE))
  cl <- call("ddspL", d = a$d, type = a$type, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("dsf.", a$type, "#", a$d, sep = "") else paste("dsf#", a$d, sep = "")
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
InitErgmTerm.dse <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","in_order"),
                      vartypes = c("numeric","character","logical"),
                      defaultvalues = list(NULL, "OTP", FALSE),
                      required = c(TRUE, FALSE,FALSE))
  cl <- call("ddspL", d = a$d, type = a$type, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("dse.", a$type, "#", a$d, sep = "") else paste("dse#", a$d, sep = "")
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
InitErgmTerm.gwdsf <- function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP",NULL,FALSE),
                      required = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))
  cl <- call("dgwdspL", decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type, alpha = a$alpha, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("gwdsf", a$type, "fixed", a$decay, sep = ".") else paste("gwdsf.fixed", a$decay, sep = ".")
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
InitErgmTerm.gwdse <- function(nw, arglist,cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP",NULL,FALSE),
                      required = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))
  cl <- call("dgwdspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type, alpha = a$alpha, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("gwdse", a$type, "fixed", a$decay, sep = ".") else paste("gwdse.fixed", a$decay, sep = ".")
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
InitErgmTerm.esf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "character,numeric","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("despL", d = a$d, type = a$type, L.base = b, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("esf(",a$base,").", a$type, "#", a$d, sep = "") else paste("esf(",a$base, ")#", a$d,sep = "")
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
InitErgmTerm.ese <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "character,numeric","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("despL",d = a$d, type = a$type, L.base = b, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("ese(",a$base,").", a$type, "#", a$d, sep = "") else paste("ese(",a$base, ")#", a$d,sep = "")
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
InitErgmTerm.gwesf <- function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "character,numeric","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL,NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("dgwespL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type, alpha = a$alpha, L.base = b, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("gwesf(",a$base,").", a$type, "fixed.", a$decay, sep = "") else paste("gwesf(",a$base, ").fixed.", a$decay,sep = "")
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
InitErgmTerm.gwese <- function(nw, arglist, cache.sp=TRUE, gw.cutoff=30,...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "character,numeric","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL,NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("dgwespL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = b, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("gwese(",a$base,").", a$type, "fixed.", a$decay, sep = "") else paste("gwese(",a$base, ").fixed.", a$decay,sep = "")
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
InitErgmTerm.nsf <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "character,numeric","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("dnspL",d = a$d, type = a$type, L.base = b, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("nsf(",a$base,").", a$type, "#", a$d, sep = "") else paste("nsf(",a$base, ")#", a$d,sep = "")
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
InitErgmTerm.nse <- function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","base", "in_order"),
                      vartypes = c("numeric","character", "character,numeric","logical"),
                      defaultvalues = list(NULL, "OTP", NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("dnspL",d = a$d, type = a$type, L.base = b, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("nse(",a$base,").", a$type, "#", a$d, sep = "") else paste("nse(",a$base, ")#", a$d,sep = "")
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
InitErgmTerm.gwnsf <- function(nw, arglist,cache.sp=TRUE, gw.cutoff=30,...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "character,numeric","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL, NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("dgwnspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = b, Ls.path= c(~`+`,~`+`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("gwnsf(",a$base,").", a$type, "fixed.", a$decay, sep = "") else paste("gwnsf(",a$base, ").fixed.", a$decay,sep = "")
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
InitErgmTerm.gwnse <- function(nw, arglist,cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay", "fixed", "cutoff","type","alpha","base","in_order"),
                      vartypes = c("numeric","logical", "numeric","character","numeric", "character,numeric","logical"),
                      defaultvalues = list(NULL,FALSE, gw.cutoff, "OTP", NULL, NULL, FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  if (is.null(a$base)) {
    b <- NULL
  } else if (a$base %in% c("+", 1)) {
    b <- ~`+`
  } else if (a$base %in% c("-", -1)) {
    b <- ~`-`
  }
  cl <- call("dgwnspL",decay = a$decay, fixed = a$fixed, cutoff = a$cutoff ,type = a$type,alpha = a$alpha, L.base = b, Ls.path= c(~`-`,~`-`), L.in_order = a$in_order)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- if(is.directed(nw)) paste("gwnse(",a$base,").", a$type, "fixed.", a$decay, sep = "") else paste("gwnse(",a$base, ").fixed.", a$decay,sep = "")
  trm
}

