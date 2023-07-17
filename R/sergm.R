#' Signed Exponential Random Graph Model (SERGM)
#'
#' The function \code{sergm} is used to fit signed exponential-family random graph models (SERGMs). The function can return a maximum pseudo-likelihood estimate, an approximate maximum likelihood estimate based on a Monte Carlo scheme, or an approximate contrastive divergenceestimate based on a similar scheme.
#'
#' @param formula An R formula object, of the form y ~ <model terms>, where y is a \code{static.sign} object. For the details on the possible <model terms>, see \link{sergm.terms}.
#' @param cons_sim Should a constraint exist that an edge can be negative and positive, default is that this is not possible.
#' @param control A list of control parameters for algorithm tuning, typically constructed with \link{control.ergm()}. Its documentation gives the the list of recognized control parameters and their meaning. The more generic utility \link{snctrl()} (StatNet ConTRoL) also provides argument completion for the available control functions and limited argument name checking.
#'
#' @return An object of class ergm that is a list consisting of coef, sample etc.
#'
#' @seealso \link{signNetwork}, \link{sergm.terms}, \link{tsergm}
#'
#' @export

sergm <- function(formula, cons_sim = T, control = control.ergm(), ...) {
  # divide formula into dependent variable (network) and independent variables
  net <- eval(formula[[2]])
  vars <- attr(terms.formula(formula),"term.labels")
  # turn signed network into multi-level network
  class(net) <- "network"
  MultiNet <- Layer(net, c(`+` = "pos",`-`= "neg"))
  # change vars to layer syntax
  vars_new <- lapply(vars, function(x) {
    shared_partners <- c("dsf", "dse", "dsm", "esf", "ese", "esm", "nsf", "nse", "nsm", "gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwesm", "gwnsf","gwnse", "gwnsm")
    degree <- c("b1degree","b2degree","degree","gwb1degree", "gwb2degree","gwdegree", "gwidegree", "gwodegree", "idegree", "odegree","isolates")
    c <- as.list(str2lang(x))
    c <- ifelse(sapply(c, is.character), paste0('"', c, '"'), c)
    temp_term <- sub("_.*", "",c[[1]])
    ### SHARED EDGES ###
      if (temp_term %in% shared_partners) {
        if (!(temp_term %in% c("gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwnsf","gwnse"))) {
          #c$d <- 1
          }
        c$L.in_order <- FALSE
        if (substr(temp_term,nchar(temp_term), nchar(temp_term)) == "f"){
          c$Ls.path <- "c(~`+`,~`+`)"
        } else if (substr(temp_term,nchar(temp_term), nchar(temp_term)) == "e"){
          c$Ls.path <- "c(~`-`,~`-`)"
        } else if (substr(temp_term,nchar(temp_term), nchar(temp_term)) == "m"){
          c$Ls.path <- "c(~`-`,~`+`)"
        } else {
          c$Ls.path <- NULL
        }
        temp_term <- gsub('e$|f$|m$', "pL", temp_term)
        if (!(as.character(temp_term) %in% c("dspL", "gwdspL"))) { #for dyadwise base can't be specified
          if (grepl("_pos", c[[1]], fixed = TRUE)) {
            c$L.base <- "~ `+`"
          } else if (grepl("_neg", c[[1]], fixed = TRUE)) {
            c$L.base <- "~ `-`"
          } else {
            c$base <- NULL
          }
        }
        c[[1]] <- temp_term
        if (length(c) > 1) {
          result <- c()
          for (i in 2:length(c)) {
            if (is.null(names(c[i])) || !nzchar(names(c[i]))) {
              result[i-1] <- c[[i]]
            } else {
              result[i-1] <- paste(names(c[i]), c[[i]], sep = "=")
            }
          }
        } else {
          result <- NULL
        }
        x <- paste(c[[1]], "(",paste(result, collapse = ", "), ")", sep="")
        }
    ### DEGREE ###
      else if (temp_term %in% degree) {
        if (temp_term == "isolates") {
          temp_term <- "degree"
          c$d <- 0
        }
          if (grepl("_pos", c[[1]], fixed = TRUE)) {
            c[[1]] <- temp_term
            if (length(c) > 1) {
              result <- c()
              for (i in 2:length(c)) {
                if (is.null(names(c[i])) || !nzchar(names(c[i]))) {
                  result[i-1] <- c[[i]]
                } else {
                  result[i-1] <- paste(names(c[i]), c[[i]], sep = "=")
                }
              }
            } else {
              result <- NULL
            }
            x <- paste("F(~",c[[1]], "(",paste(result, collapse = ", "), "), filter = ~nodefactor(~ .LayerName == '+'))", sep="")
          } else if (grepl("_neg", c[[1]], fixed = TRUE)) {
            c[[1]] <- temp_term
            if (length(c) > 1) {
              result <- c()
              for (i in 2:length(c)) {
                if (is.null(names(c[i])) || !nzchar(names(c[i]))) {
                  result[i-1] <- c[[i]]
                } else {
                  result[i-1] <- paste(names(c[i]), c[[i]], sep = "=")
                }
              }
            } else {
              result <- NULL
            }
            x <- paste("F(~",c[[1]], "(",paste(result, collapse = ", "), "), filter = ~nodefactor(~ .LayerName == '-'))", sep="")
          } else {
            c[[1]] <- temp_term
            if (length(c) > 1) {
              result <- c()
              for (i in 2:length(c)) {
                if (is.null(names(c[i])) || !nzchar(names(c[i]))) {
                  result[i-1] <- c[[i]]
                } else {
                  result[i-1] <- paste(names(c[i]), c[[i]], sep = "=")
                }
              }
            } else {
              result <- NULL
            }
            x <- paste(c[[1]], "(",paste(result, collapse = ", "), ")", sep="")
          }
      }
    ### OTHER ###
      else {
        if (grepl("_pos", c[[1]], fixed = TRUE)) {
          c[[1]] <- temp_term

          if (length(c) > 1) {
          result <- c()
          for (i in 2:length(c)) {
            if (is.null(names(c[i])) || !nzchar(names(c[i]))) {
              result[i-1] <- c[[i]]
            } else {
              result[i-1] <- paste(names(c[i]), c[[i]], sep = "=")
            }
          }
          } else {
            result <- NULL
          }

          x <- paste("L(~", c[[1]],"(", paste(result, collapse = ", ") , ")", ", ~`+`)", sep="")

        } else if (grepl("_neg", c[[1]], fixed = TRUE)) {
          c[[1]] <- temp_term

          if (length(c) > 1) {
            result <- c()
            for (i in 2:length(c)) {
              if (is.null(names(c[i])) || !nzchar(names(c[i]))) {
                result[i-1] <- c[[i]]
              } else {
                result[i-1] <- paste(names(c[i]), c[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }

          x <- paste("L(~", c[[1]],"(",paste(result, collapse = ", "), ")", ", ~`-`)", sep="")
        } else {
          x
        }
    }
  })
  # constraint that  cannot be + and - at the same time (only undirected)
  if (cons_sim == T) {
    MultiNet %ergmlhs% "constraints" <- update(MultiNet %ergmlhs% "constraints", ~ . + fixL(~`+`&`-`))
  }
  options(ergm.loglik.warn_dyads=FALSE)
  MultiFit <- ergm(as.formula(paste("MultiNet ~", paste(vars_new, collapse = "+"))),control = control, ... = ...)
  MultiFit["call"] <- format(formula)
  names(MultiFit[["MCMCtheta"]]) <- lapply(names(MultiFit[["MCMCtheta"]]), function(x) {
    a <- sub(".*~", "", x)
    if (grepl("L(pth=", x, fixed = T)) {
      b <- sub("\\..*", "", a)
      c <- sub("^[^.]+\\.", "", a)
      if (grepl("pth=(+)",x, fixed = T)) {
        b <- substring(b, 1, nchar(b)-1)
        b <- paste(b, "f", sep = "")
      } else if (grepl("pth=(-)",x, fixed = T)) {
        b <- substring(b, 1, nchar(b)-1)
        b <- paste(b, "e", sep = "")
      }
      a <- paste(b,c, sep = ".")
      if (grepl("bse=+",x, fixed = T)) {
        x <- paste(a,"_pos", sep = "")
      } else if (grepl("bse=-",x, fixed = T)) {
        x <- paste(a,"_neg", sep = "")
      } else {
        x <- a
      }
    } else {
      if (grepl("+",x, fixed = T)) {
        x <- paste(a,"_pos", sep = "")
      } else if (grepl("-",x, fixed = T)) {
        x <- paste(a,"_neg", sep = "")
      } else {
        x <- a
      }
    }
  })
  names(MultiFit[["coefficients"]]) <- lapply(names(MultiFit[["coefficients"]]), function(x) {
    a <- sub(".*~", "", x)
    if (grepl("L(pth=", x, fixed = T)) {
        b <- sub("\\..*", "", a)
        c <- sub("^[^.]+\\.", "", a)
        if (grepl("pth=(+)",x, fixed = T)) {
          b <- substring(b, 1, nchar(b)-1)
          b <- paste(b, "f", sep = "")
        } else if (grepl("pth=(-)",x, fixed = T)) {
          b <- substring(b, 1, nchar(b)-1)
          b <- paste(b, "e", sep = "")
        }
        a <- paste(b,c, sep = ".")
        if (grepl("bse=+",x, fixed = T)) {
          x <- paste(a,"_pos", sep = "")
        } else if (grepl("bse=-",x, fixed = T)) {
          x <- paste(a,"_neg", sep = "")
        } else {
          x <- a
        }
    } else {
      if (grepl("+",x, fixed = T)) {
        x <- paste(a,"_pos", sep = "")
      } else if (grepl("-",x, fixed = T)) {
        x <- paste(a,"_neg", sep = "")
      } else {
        x <- a
      }
    }
  })
  return(MultiFit)
}
