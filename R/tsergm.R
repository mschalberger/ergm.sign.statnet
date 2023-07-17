#' Temporal Signed Exponential Random Graph Model (SERGM)
#'
#' The function \code{tergm} is used for finding Temporal ERGMs’ (TERGMs) Conditional MLE (CMLE) (Krivitsky and Handcock, 2010) and Equilibrium Generalized Method of Moments Estimator (EGMME) (Krivitsky, 2009).
#'
#' @param formula An R formula object, of the form y ~ <model terms>, where y is a \code{dynamic.sign} object. For the details on the possible <model terms>, see \link{sergm.terms}.
#' @param cons_sim Should a constraint exist that an edge can be negative and positive, default is that this is not possible.
#' @param control A list of control parameters for algorithm tuning, typically constructed with \link{control.ergm()}. Its documentation gives the the list of recognized control parameters and their meaning. The more generic utility \link{snctrl()} (StatNet ConTRoL) also provides argument completion for the available control functions and limited argument name checking.
#' @param times A list of integers specifying which timepoints should be taken into account.
#'
#' @return An object of class ergm that is a list consisting of coef, sample etc.
#'
#' @seealso \link{signNetwork}, \link{sergm.terms}, \link{sergm}
#'
#' @export

tsergm <- function(formula, cons_sim = T, control = control.ergm() , times = c(0:(length(eval(formula[[2]]))-1)), ...) {
  # divide formula into dependent variable (network) and independent variables
  net <- eval(formula[[2]])
  vars <- attr(terms.formula(formula),"term.labels")
  # turn each of the networks into a Multi-level network
  NetList <- list()
  for (i in c(1:length(net))) {
    nw <- net[[i]]
    class(nw) <- "network"
    MultiNet <- Layer(nw, c(`+` = "pos",`-`= "neg"))
    NetList[[i]] <- MultiNet
  }
  # and combine them to a dynamic network
  #MultiDyn <- networkDynamic(network.list = NetList, create.TEAs = T)
  MultiDyn <- NetSeries(NetList, times = times, NA.impute = control$CMLE.NA.impute)

  # change vars to layer syntax
  vars_new <- c()
  models <- c()
  for (i in c(1:length(vars))) {

    c <- as.list(str2lang(gsub("\\+",",",gsub("~", "", vars[[i]]))))
    type <- c[[1]]
    terms <- c[-1]

    vars_new[[i]] <-  lapply(terms, function(x) {

      x <- as.list(x)
      x <- ifelse(sapply(x, is.character), paste0('"', x, '"'), x)

      shared_partners <- c("dsf", "dse", "dsm", "esf", "ese", "esm", "nsf", "nse", "nsm", "gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwesm", "gwnsf","gwnse", "gwnsm")
      degree <- c("b1degree","b2degree","degree","gwb1degree", "gwb2degree","gwdegree", "gwidegree", "gwodegree", "idegree", "odegree","isolates")

      temp_term <- sub("_.*", "", x[[1]])


      if (temp_term %in% shared_partners) {
        if (!(temp_term %in% c("gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwnsf","gwnse"))) {
          #x <- c(x,d = 1)
          }
        x <- c(x,L.in_order = FALSE)
        if (substr(temp_term,nchar(temp_term), nchar(temp_term)) == "f"){
          x <- c(x,Ls.path = "c(~`+`,~`+`)")
        } else if (substr(temp_term,nchar(temp_term), nchar(temp_term)) == "e"){
          x <- c(x,Ls.path = "c(~`-`,~`-`)")
        } else if (substr(temp_term,nchar(temp_term), nchar(temp_term)) == "m"){
          x <- c(x,Ls.path = "c(~`-`,~`+`)")
        } else {
          x <- c(x,Ls.path = NULL)
        }
        temp_term <- gsub('e$|f$|m$', "pL", temp_term)
        if (!(as.character(temp_term) %in% c("dspL", "gwdspL"))) { #for dyadwise base can't be specified
          if (grepl("_pos", x[[1]], fixed = TRUE)) {
            x <- c(x,L.base = "~ `+`")
          } else if (grepl("_neg", x[[1]], fixed = TRUE)) {
            x <- c(x,L.base = "~ `-`")
          } else {
            x <- c(x,base = NULL)
          }
        } else {
          warning("Base cannot be specified for dyadwise!")
        }
        x[[1]] <- temp_term

        if (length(x) > 1) {
          result <- c()
          for (i in 2:length(x)) {
            if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
              result[i-1] <- x[[i]]
            } else {
              result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
            }
          }
        } else {
          result <- NULL
        }

        x <- paste(x[[1]], "(",paste(result ,collapse=", " ), ")", sep="")
      }
      ### DEGREE ###
      else if (temp_term %in% degree) {
        if (temp_term == "isolates") {
          temp_term <- "degree"
          x <- c(x,d = 0)
        }
        if (grepl("_pos", x[[1]], fixed = TRUE)) {
          x[[1]] <- temp_term
          if (length(x) > 1) {
            result <- c()
            for (i in 2:length(x)) {
              if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
                result[i-1] <- x[[i]]
              } else {
                result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }
          x <- paste("F(~",x[[1]], "(",paste(result ,collapse=", " ), "), filter = ~nodefactor(~ .LayerName == '+'))", sep="")
        } else if (grepl("_neg", x[[1]], fixed = TRUE)) {
          x[[1]] <- temp_term
          if (length(x) > 1) {
            result <- c()
            for (i in 2:length(x)) {
              if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
                result[i-1] <- x[[i]]
              } else {
                result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }
          x <- paste("F(~",x[[1]], "(",paste(result ,collapse=", " ), "), filter = ~nodefactor(~ .LayerName == '-'))", sep="")
        } else {
          x[[1]] <- temp_term
          if (length(x) > 1) {
            result <- c()
            for (i in 2:length(x)) {
              if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
                result[i-1] <- x[[i]]
              } else {
                result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }
          x <- paste(x[[1]], "(",paste(result ,collapse=", " ), ")", sep="")
        }
      }
      ### OTHER ###
      else {
        if (grepl("_pos", x[[1]], fixed = TRUE)) {
          x[[1]] <- temp_term
          if (length(x) > 1) {
            result <- c()
            for (i in 2:length(x)) {
              if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
                result[i-1] <- x[[i]]
              } else {
                result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }
          x <- paste("L(~", x[[1]],"(",paste(result ,collapse=", " ), ")", ", ~`+`)", sep="")
        } else if (grepl("_neg", x[[1]], fixed = TRUE)) {
          x[[1]] <- temp_term
          if (length(x) > 1) {
            result <- c()
            for (i in 2:length(x)) {
              if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
                result[i-1] <- x[[i]]
              } else {
                result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }
          x <- paste("L(~", x[[1]],"(",paste(result ,collapse=", " ), ")", ", ~`-`)", sep="")
        } else {
          if (length(x) > 1) {
            result <- c()
            for (i in 2:length(x)) {
              if (is.null(names(x[i])) || !nzchar(names(x[i]))) {
                result[i-1] <- x[[i]]
              } else {
                result[i-1] <- paste(names(x[i]), x[[i]], sep = "=")
              }
            }
          } else {
            result <- NULL
          }
          x <- paste(x[[1]], "(" ,paste(result ,collapse=", " ), ")", sep = "")
        }
      }
    })
    models <- c(models,paste(type, "(~" ,paste(vars_new[[i]], collapse = "+"),")", sep = ""))
  }
  # estimate tergm
  if (cons_sim == T) {
    MultiDyn %ergmlhs% "constraints" <- update(MultiDyn %ergmlhs% "constraints", ~ . + fixL(~`+`&`-`))
  }
  options(ergm.loglik.warn_dyads=FALSE)
  MultiFit <- ergm(formula = as.formula(paste("MultiDyn ~", paste(models, collapse = "+"))), ... = ...#, estimate = estimate, times = times, constraints = ~. + fixL(~`+`&`-`)
                    )
  MultiFit["call"] <- format(formula)
    names(MultiFit[["MCMCtheta"]]) <- lapply(names(MultiFit[["MCMCtheta"]]), function(x) {
      a <- sub("~.*~", "~", x)
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
      a <- sub("~.*~", "~", x)
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

