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

  nw$gal$sign <- "+"

  formula <- if(is(a$formula, "formula")) list(a$formula) else a$formula
  #
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

  nw$gal$sign <- "-"

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

#' @name delese
#' @title Delayed edgewise shared enemies
#'
#' @description
#' This term adds one network statistic to the model that counts the number of
#' positive or negative edges in the current network whose endpoints had exactly
#' `d` shared enemies (common negative ties) in the previous network.
#'
#' For directed networks, different definitions of shared enemies can be specified
#' using the `type` argument. For undirected networks, only one configuration applies.
#'
#' @usage
#' # binary: delese(d = 1, base = "+", type = "OTP")
#'
#' @param d
#' Integer.
#' The exact number of shared enemies to count for edges in the current network.
#'
#' @param base
#' Character indicating which edges in the current network are used as the base:
#' `"+"` for positive ties or `"-"` for negative ties.
#'
#' @param type
#' For directed networks, the definition of shared enemies:
#' \describe{
#'   \item{"OTP"}{Outgoing two-path (\( i \to k \to j \))}
#'   \item{"ITP"}{Incoming two-path (\( j \to k \to i \))}
#'   \item{"RTP"}{Reciprocated two-path (\( i \leftrightarrow k \leftrightarrow j \))}
#'   \item{"OSP"}{Outgoing shared partner (\( i \to k, j \to k \))}
#'   \item{"ISP"}{Incoming shared partner (\( k \to i, k \to j \))}
#' }
#' Ignored for undirected networks.
#'
#' @details
#' For each edge in the current network (positive or negative, depending on `base`),
#' this term checks how many nodes were connected negatively to both endpoints
#' in the previous network.
#'
#' @concept delayed
InitErgmTerm.delese <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(
    nw, arglist,
    varnames = c("d", "base", "type"),
    vartypes = c("numeric", "character", "character"),
    defaultvalues = list(NULL, NULL, "OTP"),
    required = c(FALSE, TRUE, FALSE)
  )

  base <- match.arg(a$base, c("+", "-"))
  d <- a$d
  type <- match.arg(a$type, c("OTP", "ITP", "RTP", "OSP", "ISP"))
  is_directed <- is.directed(nw)

  # pick correct previous network
  prev_net <- nw$gal$.PrevNets[[1]]
  if (is.null(prev_net)) prev_net <- nw$gal$.PrevNet
  if (is.null(prev_net)) stop("No previous network found")

  # previous adjacency
  prev_adj <- as.matrix(get.inducedSubgraph(prev_net, which(prev_net %v% ".LayerID" == 2)))

  # shared partner matrix depending on directed type
  sp <- switch(type,
               OTP = prev_adj %*% prev_adj,
               ITP = t(prev_adj) %*% t(prev_adj),
               RTP = (prev_adj * t(prev_adj)) %*% (prev_adj * t(prev_adj)),
               OSP = t(prev_adj) %*% prev_adj,
               ISP = prev_adj %*% t(prev_adj)
  )
  diag(sp) <- 0

  # current adjacency of requested base layer
  adj <- as.matrix(get.inducedSubgraph(nw, which(nw %v% ".LayerName" == base)))
  n <- network.size(nw)
  base_nodes <- which(nw %v% ".LayerName" == base)

  # build delayed covariate matrix
  tmp <- matrix(0, nrow = n, ncol = n)
  tmp_sub <- matrix(0, length(base_nodes), length(base_nodes))
  if(!is.null(d)) tmp_sub[adj > 0 & sp == d] <- 1 else tmp_sub[adj > 0] <- sp[adj > 0]

  tmp[base_nodes, base_nodes] <- tmp_sub

  # create edgecov term
  cl <- call("edgecov", x = tmp)
  trm <- call.ErgmTerm(cl, nw, ...)

  trm$coef.names <- paste0("delese", d, base, ifelse(is_directed, paste0("~", type), ""))
  trm
}


#' @name delesf
#' @title Delayed edgewise shared friends
#'
#' @description
#' This term adds one network statistic to the model that counts the number of
#' positive or negative edges in the current network whose endpoints had exactly
#' `d` shared friends (common positive ties) in the previous network.
#'
#' For directed networks, different definitions of shared friends can be specified
#' using the `type` argument. For undirected networks, only one configuration applies.
#'
#' @usage
#' # binary: delese(d = 1, base = "+", type = "OTP")
#'
#' @param d
#' Integer.
#' The exact number of shared friends to count for edges in the current network.
#'
#' @param base
#' Character indicating which edges in the current network are used as the base:
#' `"+"` for positive ties or `"-"` for negative ties.
#'
#' @param type
#' For directed networks, the definition of shared friends:
#' \describe{
#'   \item{"OTP"}{Outgoing two-path (\( i \to k \to j \))}
#'   \item{"ITP"}{Incoming two-path (\( j \to k \to i \))}
#'   \item{"RTP"}{Reciprocated two-path (\( i \leftrightarrow k \leftrightarrow j \))}
#'   \item{"OSP"}{Outgoing shared partner (\( i \to k, j \to k \))}
#'   \item{"ISP"}{Incoming shared partner (\( k \to i, k \to j \))}
#' }
#' Ignored for undirected networks.
#'
#' @details
#' For each edge in the current network (positive or negative, depending on `base`),
#' this term checks how many nodes were connected positively to both endpoints
#' in the previous network.
#'
#' @concept delayed
InitErgmTerm.delesf <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(
    nw, arglist,
    varnames = c("d", "base", "type"),
    vartypes = c("numeric", "character", "character"),
    defaultvalues = list(NULL, NULL, "OTP"),
    required = c(FALSE, TRUE, FALSE)
  )

  base <- match.arg(a$base, c("+", "-"))
  d <- a$d
  type <- match.arg(a$type, c("OTP", "ITP", "RTP", "OSP", "ISP"))
  is_directed <- is.directed(nw)
  if (!is_directed && type != "OTP") {
    warning("Undirected network: 'type' argument ignored")
    type <- "OTP"
  }

  # pick correct previous network
  prev_net <- nw$gal$.PrevNets[[1]]
  if (is.null(prev_net)) prev_net <- nw$gal$.PrevNet
  if (is.null(prev_net)) stop("No previous network found")

  # previous adjacency
  prev_adj <- as.matrix(get.inducedSubgraph(prev_net, which(prev_net %v% ".LayerID" == 1)))

  # shared partner matrix depending on directed type
  sp <- switch(type,
               OTP = prev_adj %*% prev_adj,
               ITP = t(prev_adj) %*% t(prev_adj),
               RTP = (prev_adj * t(prev_adj) ) %*% (prev_adj * t(prev_adj)),
               OSP = t(prev_adj) %*% prev_adj,
               ISP = prev_adj %*% t(prev_adj)
  )
  diag(sp) <- 0

  # current adjacency of requested base layer
  adj <- as.matrix(get.inducedSubgraph(nw, which(nw %v% ".LayerName" == base)))
  n <- network.size(nw)
  base_nodes <- which(nw %v% ".LayerName" == base)

  # build delayed covariate matrix
  tmp <- matrix(0, nrow = n, ncol = n)
  tmp_sub <- matrix(0, length(base_nodes), length(base_nodes))
  if(!is.null(d)) tmp_sub[adj > 0 & sp == d] <- 1 else tmp_sub[adj > 0] <- sp[adj > 0]
  tmp[base_nodes, base_nodes] <- tmp_sub

  # create edgecov term
  cl <- call("edgecov", x = tmp)
  trm <- call.ErgmTerm(cl, nw, ...)

  trm$coef.names <- paste0("delesf", d, base, ifelse(is_directed, paste0("~", type), ""))
  trm
}


#' @name gwdelese
#' @title Geometrically weighted delayed edgewise shared enemies
#' @description This term calculates the number of shared enemies based on the previous network.
#' It then applies a geometric transformation to these counts to reduce the influence of large counts.
#' Specifically, if the decay parameter \code{decay} is provided, the weighting function used is
#' \deqn{f(k; decay) = 1 - (1 - e^{-decay})^k,}{f(k; d) = 1 - (1 - exp(-decay))^k,}
#' where \code{k} is the count of shared partners and \code{decay} controls how quickly the weight decreases as counts increase.
#'
#' @usage
#' gwdelese(d, base, type = "OTP")
#'
#' @param decay Numeric decay parameter controlling weighting intensity. If \code{NULL}, raw counts are used.
#' @param base Character indicating which edges in the current network are used as the base:
#' `"+"` for positive ties or `"-"` for negative ties.
#' @param type Character specifying which shared partner pattern to use; one of
#'   \code{"OTP"} (outgoing two-path),
#'   \code{"ITP"} (incoming two-path),
#'   \code{"RTP"} (reciprocated two-path),
#'   \code{"OSP"} (outgoing shared partner),
#'   \code{"ISP"} (incoming shared partner).
#'
#' @concept delayed
InitErgmTerm.gwdelese <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(
    nw, arglist,
    varnames = c("decay", "base", "type"),
    vartypes = c("numeric", "character", "character"),
    defaultvalues = list(NULL, NULL, "OTP"),
    required = c(TRUE, TRUE, FALSE)
  )

  base <- match.arg(a$base, c("+", "-"))
  decay <- a$decay
  type <- match.arg(a$type, c("OTP", "ITP", "RTP", "OSP", "ISP"))
  is_directed <- is.directed(nw)

  # pick correct previous network
  prev_net <- nw$gal$.PrevNets[[1]]
  if (is.null(prev_net)) prev_net <- nw$gal$.PrevNet
  if (is.null(prev_net)) stop("No previous network found")

  # previous adjacency matrix for selected layer nodes
  prev_adj <- as.matrix(get.inducedSubgraph(prev_net, which(prev_net %v% ".LayerID" == 2)))

  # shared partner matrix depending on directed type
  sp <- switch(type,
               OTP = prev_adj %*% prev_adj,
               ITP = t(prev_adj) %*% t(prev_adj),
               RTP = (prev_adj * t(prev_adj)) %*% (prev_adj * t(prev_adj)),
               OSP = t(prev_adj) %*% prev_adj,
               ISP = prev_adj %*% t(prev_adj)
  )
  diag(sp) <- 0

  # current adjacency of requested base layer
  adj <- as.matrix(get.inducedSubgraph(nw, which(nw %v% ".LayerName" == base)))
  n <- network.size(nw)
  base_nodes <- which(nw %v% ".LayerName" == base)

  # Build geometrically weighted covariate matrix
  tmp <- matrix(0, nrow = n, ncol = n)
  tmp_sub <- matrix(0, length(base_nodes), length(base_nodes))

  # Geometric weights: 1 - (1 - exp(-d))^sp
  # but only apply weights to edges present in adj > 0
  weights <- exp(decay)* (1 - (1 - exp(-decay))^sp)
  weights[adj <= 0] <- 0  # zero weight for non-edges

  tmp_sub <- weights

  tmp[base_nodes, base_nodes] <- tmp_sub

  # create edgecov term with weighted shared partner covariate matrix
  cl <- call("edgecov", x = tmp)
  trm <- call.ErgmTerm(cl, nw, ..., term.options = list(terms.only = TRUE))

  trm$coef.names <- paste0("gwdelese", decay, base, ifelse(is_directed, paste0("~", type), ""))
  trm
}



#' @name gwdelesf
#' @title Geometrically weighted delayed edgewise shared friends
#' @description This term calculates the number of shared friends based on the previous network.
#' It then applies a geometric transformation to these counts to reduce the influence of large counts.
#' Specifically, if the decay parameter \code{decay} is provided, the weighting function used is
#' \deqn{f(k; decay) = 1 - (1 - e^{-decay})^k,}{f(k; decay) = 1 - (1 - exp(-d))^k,}
#' where \code{k} is the count of shared partners and \code{decay} controls how quickly the weight decreases as counts increase.
#' If \code{decay} is no
#'
#' @usage
#' gwdelesf(d, base, type = "OTP")
#'
#' @param decay Numeric decay parameter controlling weighting intensity. If \code{NULL}, raw counts are used.
#' @param base Character indicating which edges in the current network are used as the base:
#' `"+"` for positive ties or `"-"` for negative ties.
#' @param type Character specifying which shared partner pattern to use; one of
#'   \code{"OTP"} (outgoing two-path),
#'   \code{"ITP"} (incoming two-path),
#'   \code{"RTP"} (reciprocated two-path),
#'   \code{"OSP"} (outgoing shared partner),
#'   \code{"ISP"} (incoming shared partner).
#'
#' @concept delayed
InitErgmTerm.gwdelesf <- function(nw, arglist, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(
    nw, arglist,
    varnames = c("decay", "base", "type"),
    vartypes = c("numeric", "character", "character"),
    defaultvalues = list(NULL, NULL, "OTP"),
    required = c(TRUE, TRUE, FALSE)
  )

  base <- match.arg(a$base, c("+", "-"))
  decay <- a$decay
  type <- match.arg(a$type, c("OTP", "ITP", "RTP", "OSP", "ISP"))
  is_directed <- is.directed(nw)

  # pick correct previous network
  prev_net <- nw$gal$.PrevNets[[1]]
  if (is.null(prev_net)) prev_net <- nw$gal$.PrevNet
  if (is.null(prev_net)) stop("No previous network found")

  # previous adjacency matrix for selected layer nodes
  prev_adj <- as.matrix(get.inducedSubgraph(prev_net, which(prev_net %v% ".LayerID" == 1)))

  # shared partner matrix depending on directed type
  sp <- switch(type,
               OTP = prev_adj %*% prev_adj,
               ITP = t(prev_adj) %*% t(prev_adj),
               RTP = (prev_adj * t(prev_adj)) %*% (prev_adj * t(prev_adj)),
               OSP = t(prev_adj) %*% prev_adj,
               ISP = prev_adj %*% t(prev_adj)
  )
  diag(sp) <- 0

  # current adjacency of requested base layer
  adj <- as.matrix(get.inducedSubgraph(nw, which(nw %v% ".LayerName" == base)))
  n <- network.size(nw)
  base_nodes <- which(nw %v% ".LayerName" == base)

  # Build geometrically weighted covariate matrix
  tmp <- matrix(0, nrow = n, ncol = n)
  tmp_sub <- matrix(0, length(base_nodes), length(base_nodes))

  # Geometric weights: 1 - (1 - exp(-d))^sp
  # but only apply weights to edges present in adj > 0
  weights <-exp(decay)* (1 - (1 - exp(-decay))^sp)
  weights[adj <= 0] <- 0  # zero weight for non-edges

  tmp_sub <- weights

  tmp[base_nodes, base_nodes] <- tmp_sub

  # create edgecov term with weighted shared partner covariate matrix
  cl <- call("edgecov", x = tmp)
  trm <- call.ErgmTerm(cl, nw, ..., terms.only = TRUE)

  trm$coef.names <- paste0("gwdelesf", decay, base, ifelse(is_directed, paste0("~", type), ""))
  trm
}


#' @templateVar name delrecip
#' @title Delayed reciprocity
#' @description For the current network layer `base`, this term equals 1 for each directed edge
#' i->j currently present where the reverse edge j->i was present in the previous network's same
#' layer. The previous network used is the one indexed by the current `GroupID` (i.e., the behaviour
#' is the same as the prior `lag=1` implementation). The term is provided as an edgecov (1/0).
#'
#' @usage
#' # binary: delrecip(base)
#'
#' @param base character or numeric name/identifier of the layer to examine in the current network.
#'
#' @concept delayed
InitErgmTerm.delrecip <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = "base",
                      vartypes = "character",
                      defaultvalues = list(NULL),
                      required = FALSE)

  base <- a$base
  nw1 <- nw

  prev_net <- nw$gal$.PrevNets[[1]]
  if (is.null(prev_net)) prev_net <- nw$gal$.PrevNet
  if (is.null(prev_net)) stop("No previous network found")

  # adjacency of previous net restricted to requested layer
  if (!is.null(base)) {
  prev_net <- get.inducedSubgraph(prev_net, which(prev_net %v% ".LayerName" == base))
  nw1 <- get.inducedSubgraph(nw, which(nw %v% ".LayerName" == base))
  }
  prev_adj <- as.matrix(prev_net)
  curr_adj <- as.matrix(nw1)

  # reciprocity: 1 if current i->j and previous j->i
  result_matrix <- (curr_adj == 1 & t(prev_adj) == 1) * 1

  n <- network.size(nw)
  base_nodes <- which(nw %v% ".LayerName" == base)

  tmp <- matrix(0, nrow = n, ncol = n)
  tmp[base_nodes, base_nodes] <- result_matrix

  cl <- call("edgecov", x = tmp)
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- paste("delrecip")
  trm
}


#' @templateVar name delnodematch
#' @title Delayed node matching on attribute (lag-1)
#' @description Create a nodematch term where node attributes come from the previous network's
#' node attribute `attr`. The previous network used is the one indexed by the current `GroupID`
#' (equivalent to `lag=1` previously). This constructs a temporary node attribute `delnodecov_<attr>`
#' on the current network (copied from the previous net) and calls `nodematch` on that attribute.
#'
#' @usage
#' # binary: delnodematch(attr)
#'
#' @param attr character attribute name to copy from the previous network into the current.
#'
#' @concept delayed
InitErgmTerm.delnodematch <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  prev_net <- nw$gal$.PrevNets[[1]]
  if (is.null(prev_net)) prev_net <- nw$gal$.PrevNet
  if (is.null(prev_net)) stop("No previous network found")

  # if (!is.null(nw$gal$sign)) {
  #   prev_net <- get.inducedSubgraph(prev_net, which(prev_net %v% ".LayerName" == nw$gal$sign
  #   ))
  #   }

  # copy attribute values from previous net to a temporary attribute on current net
  nw %v% paste("delnodecov", a$attr, sep = "_") <- prev_net %v% a$attr

  # Call nodematch on the temporary attribute
  cl <- call("nodematch", paste("delnodecov", a$attr, sep = "_"))
  trm <- call.ErgmTerm(cl, nw, ...)
  trm$coef.names <- paste("delnodematch", a$attr, sep = "_")
  trm
}


#InitErgmTerm..layer.net <- ergm.multi:::InitErgmTerm..layer.net
#body(InitErgmTerm..layer.net)[[7]] <- quote(list(name="_layer_net", coef.names=c(), iinputs=c(unlist(.block_vertexmap(nw, ".LayerID", FALSE)), if(is.directed(nw)) sapply(nwl, function(nw) (nw%v% ".undirected")[1]), ll), dependence=dependence))
#assignInNamespace("InitErgmTerm..layer.net", InitErgmTerm..layer.net, ns="ergm.multi")

# uncombine_network <- ergm.multi:::uncombine_network
# body(uncombine_network)[[6]][[4]][[4]][[4]][[2]] <- quote({
#   if ("combined_networks" %in% class(nwl[[i]])) {
#     class(nwl[[i]]) <- class(nwl[[i]])[-seq_len(min(which(class(nwl[[i]]) == "combined_networks")))]
#   }
# })
# assignInNamespace("uncombine_network", uncombine_network, ns="ergm.multi")
#
# InitErgmTerm.N <- ergm.multi:::InitErgmTerm.N
# body(InitErgmTerm.N)[[7]] <- quote(nattrs <- as_tibble.combined_networks(nw, unit = "networks", .NetworkID = a$.NetworkID,
#                                                         .NetworkName = a$.NetworkName))
# assignInNamespace("InitErgmTerm.N", InitErgmTerm.N, ns="ergm.multi")
#
# patched_uncombine <-
#   function(nw,
#            split.vattr = nw %n% ".blockID.vattr",
#            names.vattr = nw %n% ".blockName.vattr",
#            use.subnet.cache = FALSE) {
#     if (!is(nw, "combined_networks"))
#       stop("Specified network was not constructed by combine_networks().")
#
#     tmp <- ergm.multi:::.pop_vattrv(nw, split.vattr)
#     nw <- tmp$nw
#     f <- tmp$vattr
#
#     if (!is.null(names.vattr)) {
#       tmp <- ergm.multi:::.pop_vattrv(nw, names.vattr)
#       nw <- tmp$nw
#       nwnames <- tmp$vattr
#     }
#
#     # if (use.subnet.cache && ".subnetcache" %in% list.network.attributes(nw) &&
#     #     names(nw %n% ".subnetcache") == split.vattr) {
#     #   nwl <- (nw %n% ".subnetcache")[[split.vattr]]
#     # } else {
#       # just split, skip the for-loop
#       nwl <- split(nw, f)
#     #}
#
#     if (!is.null(names.vattr))
#       names(nwl) <- unique(nwnames)
#
#     class(nwl) <- c("network.list", class(nwl))
#     nwl
#   }
#

