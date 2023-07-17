#' Draw from the Distribution of a Temporal Signed Exponential Family Random Graph Model
#'
#' \code{sim_tsergm} is used to draw from temporal signed exponential family random network models. See \link{tsergm} for more information on these models. The method for tsergm objects inherits the model, the coefficients, the response attribute, the reference, the constraints, and most simulation parameters from the model fit, unless overridden by passing them explicitly. Unless overridden, the simulation is initialized with either a random draw from near the fitted model saved by \code{tsergm()}.
#'
#'
#' @param object Either a formula or a fitted tsergm. The formula should be of the form y ~ <model terms>, where y is a signed network of the class \code{dynamic.sign}.
#' @param nsim Number of networks to be randomly drawn from the given distribution on the set of all networks, returned by the Metropolis-Hastings algorithm.
#' @param seed Seed value (integer) for the random number generator. See \link{set.seed}.
#' @param coef Vector of parameter values for the model from which the sample is to be drawn.
#' @param ... Further arguments passed to or used by methods.
#'
#' @return A dynamic signed network or a list of dynamic signed networks (if nsim > 1) of class \code{dynamic.sign}.
#'
#' @seealso \link{signNetwork}, \link{tsergm}, \link{sergm.terms}, \link{sim_tsergm}
#'
#' @export


sim_tsergm <- function(object, nsim = 1, seed = NULL, coef = NULL, ...) {

  if (class(object) == "formula") {
    net <- eval(formula[[2]])
    vars <- attr(terms.formula(formula),"term.labels")
    # turn each of the networks into a Multi-level network
    if (class(net) == "dynamic.sign"||class(net) == "network.list" || class(net) == "list") {
    NetList <- list()
    for (i in c(1:length(net))) {
      nw <- net[[i]]
      class(nw) <- "network"
      MultiNet <- Layer(nw, c(`+` = "pos",`-`= "neg"))
      NetList[[i]] <- MultiNet
    }
    # and combine them to a dynamic network
    #MultiDyn <- networkDynamic(network.list = NetList, create.TEAs = T)
    MultiDyn <- NetSeries(NetList, NA.impute = control.ergm()$CMLE.NA.impute)
    } else {
      class(net) <- "network"
      MultiDyn <- Layer(net, c(`+` = "pos",`-`= "neg"))
    }
    # change vars to layer syntax
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
          temp_term <- gsub('e$|f$', "pL", temp_term)
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
    MultiDyn %ergmlhs% "constraints" <- update(MultiDyn %ergmlhs% "constraints", ~ . + fixL(~`+`&`-`))
    sim <- simulate(as.formula(paste("MultiDyn ~", paste(models, collapse = "+"))), nsim = nsim, seed = seed, coef = coef, dynamic = TRUE)
    if (nsim > 1) {
      result <- list()
      for(i in c(1:length(sim))) {
        nw <- sim[[i]]
        sim_matrix <- as.matrix.network(nw)
        netid <- get.vertex.attribute(nw, attrname = ".NetworkID")

        res <- lapply(unique(netid), function(i) {
          end_index <- sum(netid <= i)
          start_index <- sum(netid <= i-1) + 1
          net <- sim_matrix[start_index:end_index, start_index:end_index]
          n <- ncol(net)
          net[1:(n/2),1:(n/2)] + net[(n/2+1):n,(n/2+1):n]*-1
        })
        result[[i]] <- signNetwork(res, matrix.type = "adjacency", ... = ...)
      }
    } else {
      sim_matrix <- as.matrix.network(sim)
      netid <- get.vertex.attribute(sim, attrname = ".NetworkID")

      res <- lapply(unique(netid), function(i) {
        end_index <- sum(netid <= i)
        start_index <- sum(netid <= i-1) + 1
        net <- sim_matrix[start_index:end_index, start_index:end_index]
        n <- ncol(net)
        net[1:(n/2),1:(n/2)] + net[(n/2+1):n,(n/2+1):n]*-1
      })
      result <- signNetwork(res, matrix.type = "adjacency", ... = ...)
    }
    return(result)

  } else if (class(object) == "ergm") {
    if (is.null(coef)) {
      coef <- coefficients(object)
    }

  sim <- simulate(object, nsim = nsim, seed = seed, coef = coef)

  if (nsim > 1) {
    result <- list()
    for(i in c(1:length(sim))) {
      nw <- sim[[i]]
      sim_matrix <- as.matrix.network(nw)
      netid <- get.vertex.attribute(nw, attrname = ".NetworkID")

      res <- lapply(unique(netid), function(i) {
        end_index <- sum(netid <= i)
        start_index <- sum(netid <= i-1) + 1
        net <- sim_matrix[start_index:end_index, start_index:end_index]
        n <- ncol(net)
        net[1:(n/2),1:(n/2)] + net[(n/2+1):n,(n/2+1):n]*-1
      })
      result[[i]] <- signNetwork(res, matrix.type = "adjacency", ... = ...)
    }
  } else {
    sim_matrix <- as.matrix.network(sim)
    netid <- get.vertex.attribute(sim, attrname = ".NetworkID")

    res <- lapply(unique(netid), function(i) {
      end_index <- sum(netid <= i)
      start_index <- sum(netid <= i-1) + 1
      net <- sim_matrix[start_index:end_index, start_index:end_index]
      n <- ncol(net)
      net[1:(n/2),1:(n/2)] + net[(n/2+1):n,(n/2+1):n]*-1
    })
    result <- signNetwork(res, matrix.type = "adjacency", ... = ...)
  }
  return(result)
  }
}

