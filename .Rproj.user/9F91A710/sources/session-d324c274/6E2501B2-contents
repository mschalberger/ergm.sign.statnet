#' Calculation of network or graph statistics or other attributes specified in a formula
#'
#' This function computes summaries of the object on the left-hand side (LHS) of the formula that are specified by its right-hand side (RHS). If a network is given as the LHS and \link{sergm.terms} is on the RHS, it computes the sufficient statistics associated with those terms.
#'
#' @param formula A formula having as its LHS a \code{static.sign} or \code{dynamic.sign} object to be summarized using a formula.
#'
#' @seealso \link{signnet}, \link{sergm.terms}
#'
#' @return A vector of statistics specified in the RHS of the formula.
#'
#' @export

count <- function(formula) {
  # divide formula into dependent variable (network) and independent variables
  net <- eval(formula[[2]])
  vars <- attr(terms.formula(formula),"term.labels")

  if (class(net) == "dynamic.sign" & (grepl("Cross", vars) || grepl("Change", vars) ||
      grepl("Form", vars) || grepl("Diss", vars))) {

    MultiNet <- list()
    for (i in c(1:length(net))) {
      nw <- net[[i]]
      class(nw) <- "network"
      Multi <- Layer(nw, c(`+` = "pos",`-`= "neg"))
      MultiNet[[i]] <- Multi
    }
    # and combine them to a dynamic network
    MultiNet <- NetSeries(MultiNet, NA.impute = control.ergm()$CMLE.NA.impute)
    #class(MultiNet) <- "network.list"

    vars_new <- c()
    models <- c()
    for (i in c(1:length(vars))) {

      c <- as.list(str2lang(gsub("\\+",",",gsub("~", "", vars[[i]]))))
      type <- c[[1]]
      terms <- c[-1]

      vars_new[[i]] <-  lapply(terms, function(x) {

        x <- as.list(x)

        shared_partners <- c("dsf", "dse", "dsm", "esf", "ese", "nsf", "nse", "gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwnsf","gwnse")
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
          }
          x[[1]] <- temp_term
          x <- paste(x[[1]], "(",paste(names(x[-1]),x[-1],sep="=",collapse=", " ), ")", sep="")
        }
        ### DEGREE ###
        else if (temp_term %in% degree) {
          if (temp_term == "isolates") {
            temp_term <- "degree"
            x <- c(x,d = 0)
          }
          if (grepl("_pos", x[[1]], fixed = TRUE)) {
            x[[1]] <- temp_term
            x <- paste("F(~",x[[1]], "(",paste(names(x[-1]),x[-1],sep="=",collapse=", " ), "), filter = ~nodefactor(~ .LayerName == '+'))", sep="")
          } else if (grepl("_neg", x[[1]], fixed = TRUE)) {
            x[[1]] <- temp_term
            x <- paste("F(~",x[[1]], "(",paste(names(x[-1]),x[-1],sep="=",collapse=", " ), "), filter = ~nodefactor(~ .LayerName == '-'))", sep="")
          } else {
            x[[1]] <- temp_term
            x <- paste(x[[1]], "(",paste(names(x[-1]),x[-1],sep="=",collapse=", " ), ")", sep="")
          }
        }
        ### OTHER ###
        else {
          if (grepl("_pos", x[[1]], fixed = TRUE)) {
            x[[1]] <- temp_term
            x <- paste("L(~", x[[1]],"(",paste(names(x[-1]),x[-1],sep="=",collapse=", " ), ")", ", ~`+`)", sep="")
          } else if (grepl("_neg", x[[1]], fixed = TRUE)) {
            x[[1]] <- temp_term
            x <- paste("L(~", x[[1]],"(",paste(names(x[-1]),x[-1],sep="=",collapse=", " ), ")", ", ~`-`)", sep="")
          } else {
            x
          }
        }
      })
      models <- c(models,paste(type, "(~" ,paste(vars_new[[i]], collapse = "+"),")", sep = ""))
    }
    MultiFit <- summary_formula(object = as.formula(paste("MultiNet ~", paste(models, collapse = "+"))))
    names(MultiFit) <- lapply(names(MultiFit), function(x) {
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

  }
  else {

  if(class(net) == "static.sign" | class(net) == "network") {
  # turn signed network into multi-level network
  class(net) <- "network"
  MultiNet <- Layer(net, c(`+` = "pos",`-`= "neg"))
  }
  else if (class(net) == "dynamic.sign") {
      MultiNet <- list()
    for (i in c(1:length(net))) {
      nw <- net[[i]]
      class(nw) <- "network"
      Multi <- Layer(nw, c(`+` = "pos",`-`= "neg"))
      MultiNet[[i]] <- Multi
    }
      class(MultiNet) <- "network.list"
    }
    # and combine them to a dynamic network
    #MultiNet <- NetSeries(NetList, NA.impute = control.ergm()$CMLE.NA.impute)


  # change vars to layer syntax
  vars_new <- lapply(vars, function(x) {
    shared_partners <- c("dsf", "dse", "dsm", "esf", "ese", "nsf", "nse", "gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwnsf","gwnse")
    degree <- c("b1degree","b2degree","degree","gwb1degree", "gwb2degree","gwdegree", "gwidegree", "gwodegree", "idegree", "odegree","isolates")
    c <- as.list(str2lang(x))
    temp_term <- sub("_.*", "",c[[1]])
    ### SHARED EDGES ###
    if (temp_term %in% shared_partners) {
      if (!(temp_term %in% c("gwdsf", "gwdse", "gwdsm","gwesf","gwese", "gwnsf","gwnse"))) {
        c$d <- 1}
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
      x <- paste(c[[1]], "(",paste(names(c[-1]),c[-1],sep="=",collapse=", " ), ")", sep="")
    }
    ### DEGREE ###
    else if (temp_term %in% degree) {
      if (temp_term == "isolates") {
        temp_term <- "degree"
        c$d <- 0
      }
      if (grepl("_pos", c[[1]], fixed = TRUE)) {
        c[[1]] <- temp_term
        x <- paste("F(~",c[[1]], "(",paste(names(c[-1]),c[-1],sep="=",collapse=", " ), "), filter = ~nodefactor(~ .LayerName == '+'))", sep="")
      } else if (grepl("_neg", c[[1]], fixed = TRUE)) {
        c[[1]] <- temp_term
        x <- paste("F(~",c[[1]], "(",paste(names(c[-1]),c[-1],sep="=",collapse=", " ), "), filter = ~nodefactor(~ .LayerName == '-'))", sep="")
      } else {
        c[[1]] <- temp_term
        x <- paste(c[[1]], "(",paste(names(c[-1]),c[-1],sep="=",collapse=", " ), ")", sep="")
      }
    }
    ### OTHER ###
    else {
      if (grepl("_pos", c[[1]], fixed = TRUE)) {
        c[[1]] <- temp_term
        x <- paste("L(~", c[[1]],"(",paste(names(c[-1]),c[-1],sep="=",collapse=", " ), ")", ", ~`+`)", sep="")
      } else if (grepl("_neg", c[[1]], fixed = TRUE)) {
        c[[1]] <- temp_term
        x <- paste("L(~", c[[1]],"(",paste(names(c[-1]),c[-1],sep="=",collapse=", " ), ")", ", ~`-`)", sep="")
      } else {
        x
      }
    }
  })

  MultiFit <- summary_formula(as.formula(paste("MultiNet ~", paste(vars_new, collapse = "+"))))
  if(class(net) == "static.sign" | class(net) == "network") {
  names(MultiFit) <- lapply(names(MultiFit), function(x) {
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
  } else if (class(net) == "dynamic.sign") {
    colnames(MultiFit) <- lapply(colnames(MultiFit), function(x) {
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
    rownames(MultiFit) <- c(1:length(net))
  }
  }
  return(MultiFit)
}

