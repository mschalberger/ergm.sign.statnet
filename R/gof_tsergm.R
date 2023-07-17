#' Conduct Goodness-of-Fit Diagnostics on a Temporal Signed Exponential Family Random Graph Model
#'
#' This function calculates the goodness-of-fit (GoF) for a fitted TSERGM model by simulating new networks using the same model and comparing the average observed statistics of the original network over the timepoints with those of the simulated networks. The statistics used in the comparison are the positive and negative degree distributions, edge-wise shared enemies distributions for positive and negative edges, and edge-wise shared friend distributions for positive and negative edges.
#'
#' @param model A fitted TSERGM model.
#' @param nsim An integer representing the number of simulated networks to generate. Defaults to 200.
#'
#' @return Plots 6 diagnostics for the goodness-of-fit of temporal signed exponential family random graph models.
#'
#' @seealso \link{tsergm}, \link{gof_sergm}
#'
#' @export

gof_tsergm <- function(model, nsim = 200) {

  suppressWarnings({
  # simulate tsergms
  sim <- sim_tsergm(model, nsim = nsim)

  sim_matrix <- as.matrix.network(model[["network"]])
  netid <- get.vertex.attribute(model[["network"]], attrname = ".NetworkID")

  res <- lapply(unique(netid), function(i) {
    end_index <- sum(netid <= i)
    start_index <- sum(netid <= i-1) + 1
    net <- sim_matrix[start_index:end_index, start_index:end_index]
    n <- ncol(net)
    net[1:(n/2),1:(n/2)] + net[(n/2+1):n,(n/2+1):n]*-1
  })

  net_list <- signNetwork(res, matrix.type = "adjacency")
  class(net_list) <- "network"

  n <- e <- degree_pos <- degree_neg <- ese_neg <- ese_pos <- esf_neg <- esf_pos <- c()
  for (i in c(1:length(net_list))) {
    nw <- net_list[[i]]
    n[[i]] <- network.size(nw)
    MultiNet <- Layer(nw, c(`+` = "pos",`-`= "neg"))

    e[[i]] <- summary_formula(MultiNet ~ edges)

    degree_pos[[i]] <- summary_formula(MultiNet ~ L(~ degree(c(0:(n[[i]]-1))), ~ `+`)) #/ n[[i]]
    names(degree_pos[[i]]) <- c(0:(n[[i]]-1))

    degree_neg[[i]] <- summary_formula(MultiNet ~ L(~ degree(c(0:(n[[i]]-1))), ~ `-`)) #/ n[[i]]
    names(degree_neg[[i]]) <- c(0:(n[[i]]-1))

    ese_neg[[i]] <- summary_formula(MultiNet ~ espL(c(0:(n[[i]]-2)), L.base=~`-`, Ls.path=c(~`-`,~`-`))) #/ e[[i]]
    names(ese_neg[[i]]) <- c(0:(n[[i]]-2))

    ese_pos[[i]] <- summary_formula(MultiNet ~ espL(c(0:(n[[i]]-2)), L.base=~`+`, Ls.path=c(~`-`,~`-`))) #/ e[[i]]
    names(ese_pos[[i]]) <- c(0:(n[[i]]-2))

    esf_neg[[i]] <- summary_formula(MultiNet ~ espL(c(0:(n[[i]]-2)), L.base=~`-`, Ls.path=c(~`+`,~`+`))) #/ e[[i]]
    names(esf_neg[[i]]) <- c(0:(n[[i]]-2))

    esf_pos[[i]] <- summary_formula(MultiNet ~ espL(c(0:(n[[i]]-2)), L.base=~`+`, Ls.path=c(~`+`,~`+`))) #/ e[[i]]
    names(esf_pos[[i]]) <- c(0:(n[[i]]-2))
  }

  n <- Reduce("+", n)/length(net_list)
  e <- Reduce("+", e)/length(net_list)
  degree_pos <- Reduce("+", degree_pos)/length(net_list)
  degree_neg <- Reduce("+", degree_neg)/length(net_list)
  ese_neg <- Reduce("+", ese_neg)/length(net_list)
  ese_pos <- Reduce("+", ese_pos)/length(net_list)
  esf_neg <- Reduce("+", esf_neg)/length(net_list)
  esf_pos <- Reduce("+", esf_pos)/length(net_list)

  sim_degree_pos <- sim_degree_neg <- net_size <- edges <- sim_ese_neg <- sim_ese_pos <- sim_esf_neg <- sim_esf_pos <- list()
  for (i in c(1:length(sim))) {

    tmp_net <- sim[[i]]

    sim_degree_pos[[i]] <- list()
    sim_degree_neg[[i]] <- list()
    net_size[[i]] <- list()
    edges[[i]] <- list()
    sim_ese_neg[[i]] <- list()
    sim_ese_pos[[i]] <- list()
    sim_esf_neg[[i]] <- list()
    sim_esf_pos[[i]] <- list()

    for (j in c(1:length(tmp_net))) {
      net_size[[i]][[j]] <- network.size(tmp_net[[j]])
      tmp <- Layer(tmp_net[[j]], c(`+` = "pos",`-`= "neg"))
      edges[[i]][[j]] <- summary_formula(tmp ~ edges)
      sim_degree_pos[[i]][[j]] <- summary_formula(tmp ~ L(~ degree(c(0:(net_size[[i]][[j]]-1))), ~ `+`)) #/ net_size[[i]][[j]]
      sim_degree_neg[[i]][[j]] <- summary_formula(tmp ~ L(~ degree(c(0:(net_size[[i]][[j]]-1))), ~ `-`)) #/ net_size[[i]][[j]]
      sim_ese_neg[[i]][[j]] <- summary_formula(tmp ~ espL(c(0:(net_size[[i]][[j]]-2)), L.base=~`-`, Ls.path=c(~`-`,~`-`))) #/ edges[[i]][[j]]
      sim_ese_pos[[i]][[j]] <- summary_formula(tmp ~ espL(c(0:(net_size[[i]][[j]]-2)), L.base=~`+`, Ls.path=c(~`-`,~`-`))) #/ edges[[i]][[j]]
      sim_esf_neg[[i]][[j]] <- summary_formula(tmp ~ espL(c(0:(net_size[[i]][[j]]-2)), L.base=~`-`, Ls.path=c(~`+`,~`+`))) #/ edges[[i]][[j]]
      sim_esf_pos[[i]][[j]] <-  summary_formula(tmp ~ espL(c(0:(net_size[[i]][[j]]-2)), L.base=~`+`, Ls.path=c(~`+`,~`+`))) #/ edges[[i]][[j]]
    }
  }

  list_dp <- lapply(sim_degree_pos, function (x) {
    x <- Reduce("+", x)/length(net_list)
  })
  list_dp <- bind_rows(list_dp)
  colnames(list_dp) <- seq(0, ncol(list_dp))
  list_dp <- as.data.frame(list_dp)

  list_dn <- lapply(sim_degree_neg, function (x) {
    x <- Reduce("+", x)/length(net_list)
  })
  list_dn <- bind_rows(list_dn)
  colnames(list_dn) <- seq(0, ncol(list_dn))
  list_dn <- as.data.frame(list_dn)

  list_sen <- lapply(sim_ese_neg, function (x) {
    x <- Reduce("+", x)/length(net_list)
  })
  list_sen <- bind_rows(list_sen)
  colnames(list_sen) <- seq(0, ncol(list_sen))
  list_sen <- as.data.frame(list_sen)

  list_sep <- lapply(sim_ese_pos, function (x) {
    x <- Reduce("+", x)/length(net_list)
  })
  list_sep <- bind_rows(list_sep)
  colnames(list_sep) <- seq(0, ncol(list_sep))
  list_sep <- as.data.frame(list_sep)

  list_sfn <- lapply(sim_esf_neg, function (x) {
    x <- Reduce("+", x)/length(net_list)
  })
  list_sfn <- bind_rows(list_sfn)
  colnames(list_sfn) <- seq(0, ncol(list_sfn))
  list_sfn <- as.data.frame(list_sfn)

  list_sfp <- lapply(sim_esf_pos, function (x) {
    x <- Reduce("+", x)/length(net_list)
  })
  list_sfp <- bind_rows(list_sfp)
  colnames(list_sfp) <- seq(0, ncol(list_sfp))
  list_sfp <- as.data.frame(list_sfp)

  # Extract max index
  get_max_indices <- function(x) {
    # Get column indices where the absolute sum is not zero
    col_indices <- which(sapply(x, function(y) sum(abs(y))) != 0)
    # Get the index of the maximum value
    max_indices <- max(col_indices)
    # Return the max indices
    return(max_indices)
  }


  boxplot(list_dp,
          names = names(list_dp),
          xlab = "Positive Degrees",
          ylab = "Average number of observations",
          ylim = c(0,max(unlist(list_dp), unlist(degree_pos), na.rm = T)),
          xlim = c(0, max(get_max_indices(list_dp), get_max_indices(degree_pos))+ 1))
  lines(degree_pos, lwd = 3)

  boxplot(list_dn,
          names = names(list_dn),
          xlab = "Negative degrees",
          ylab = "Average number of observations",
          ylim = c(0,max(unlist(list_dn), unlist(degree_neg), na.rm = T)),
          xlim = c(0, max(get_max_indices(list_dn), get_max_indices(degree_neg))+ 1))
  lines(degree_neg, lwd = 3)

  boxplot(list_sen,
          names = names(list_sen),
          xlab = "k−Edgewise-Shared Enemies (-)",
          ylab = "Average number of observations",
          ylim = c(0,max(unlist(list_sen), unlist(ese_neg), na.rm = T)),
          xlim = c(0, max(get_max_indices(list_sen), get_max_indices(ese_neg))+ 1))
  lines(ese_neg, lwd = 3)

  boxplot(list_sep,
          names = names(list_sep),
          xlab = "k−Edgewise-Shared Enemies (+)",
          ylab = "Average number of observations",
          ylim = c(0,max(unlist(list_sep), unlist(ese_pos), na.rm = T)),
          xlim = c(0, max(get_max_indices(list_sep), get_max_indices(ese_pos))+ 1))
  lines(ese_pos, lwd = 3)

  boxplot(list_sfn,
          names = names(list_sfn),
          xlab = "k−Edgewise-Shared Friends (-)",
          ylab = "Average number of observations",
          ylim = c(0,max(unlist(list_sfn), unlist(esf_neg), na.rm = T)),
          xlim = c(0, max(get_max_indices(list_sfn), get_max_indices(esf_neg))+ 1))
  lines(esf_neg, lwd = 3)

  boxplot(list_sfp,
          names = names(list_sfp),
          xlab = "k−Edgewise-Shared Friends (+)",
          ylab = "Average number of observations",
          ylim = c(0,max(unlist(list_sfp), unlist(esf_pos), na.rm = T)),
          xlim = c(0, max(get_max_indices(list_sfp), get_max_indices(esf_pos)) + 1))
  lines(esf_pos, lwd = 3)

  })
}
