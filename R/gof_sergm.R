#' Conduct Goodness-of-Fit Diagnostics on a Signed Exponential Family Random Graph Model
#'
#' This function calculates the goodness-of-fit (GoF) for a fitted SERGM model by simulating new networks using the same model and comparing the statistics of the original network with those of the simulated networks. The statistics used in the comparison are the positive and negative degree distributions, edge-wise shared enemies distributions for positive and negative edges, and edge-wise shared friend distributions for positive and negative edges.
#'
#' @param model A fitted SERGM model.
#' @param nsim An integer representing the number of simulated networks to generate. Defaults to 200.
#'
#' @return Plots 6 diagnostics for the goodness-of-fit of signed exponential family random graph models.
#'
#' @seealso \link{sergm}, \link{gof_tsergm}
#'
#' @export

gof_sergm <- function(model, nsim = 200) {

  # simulate sergms
  sim <- sim_sergm(model, nsim = nsim)

  # extract network from model
  matrix <- as.matrix.network(model[["network"]])


  pos <- as.sociomatrix(model[["network"]][["gal"]][[".subnetcache"]][[".LayerID"]][["+"]])
  neg <- as.sociomatrix(model[["network"]][["gal"]][[".subnetcache"]][[".LayerID"]][["-"]])*-1
  comb <- pos + neg
  net <- signNetwork(mat = comb, matrix.type = "adjacency")
  class(net) <- "network"

  n <- network.size(net)
  MultiNet <- Layer(net, c(`+` = "pos",`-`= "neg"))
  e <- summary_formula(MultiNet ~ edges)

  # create network statistics
  degree_pos <- summary_formula(MultiNet ~ L(~ degree(c(0:(n-1))), ~ `+`)) / n
  names(degree_pos) <- c(0:(n-1))

  degree_neg <- summary_formula(MultiNet ~ L(~ degree(c(0:(n-1))), ~ `-`)) / n
  names(degree_neg) <- c(0:(n-1))

  ese_neg <- summary_formula(MultiNet ~ espL(c(0:(n-2)), L.base=~`-`, Ls.path=c(~`-`,~`-`))) / e
  names(ese_neg) <- c(0:(n-2))

  ese_pos <- summary_formula(MultiNet ~ espL(c(0:(n-2)), L.base=~`+`, Ls.path=c(~`-`,~`-`))) / e
  names(ese_pos) <- c(0:(n-2))

  esf_neg <- summary_formula(MultiNet ~ espL(c(0:(n-2)), L.base=~`-`, Ls.path=c(~`+`,~`+`))) / e
  names(esf_neg) <- c(0:(n-2))

  esf_pos <- summary_formula(MultiNet ~ espL(c(0:(n-2)), L.base=~`+`, Ls.path=c(~`+`,~`+`))) / e
  names(esf_pos) <- c(0:(n-2))

  # create network statistics for simulated networks
   sim_degree_pos <- sim_degree_neg <- net_size <- edges <- sim_ese_neg <- sim_ese_pos <- sim_esf_neg <- sim_esf_pos <- c()
   for (i in c(1:nsim)) {
     class(sim[[i]]) <- "network"
     tmp <- Layer(sim[[i]], c(`+` = "pos",`-`= "neg"))
     net_size[i] <- network.size(sim[[i]])
     edges[i] <- summary_formula(tmp ~ edges)
     sim_degree_pos[[i]] <- summary_formula(tmp ~ L(~ degree(c(0:(network.size(sim[[i]])-1))), ~ `+`))/net_size[i]
     sim_degree_neg[[i]] <- summary_formula(tmp ~ L(~ degree(c(0:(network.size(sim[[i]])-1))), ~ `-`))/net_size[i]
     sim_ese_neg[[i]] <- summary_formula(tmp ~ espL(c(0:(network.size(sim[[i]])-2)), L.base=~`-`, Ls.path=c(~`-`,~`-`)))/edges[i]
     sim_ese_pos[[i]] <- summary_formula(tmp ~ espL(c(0:(network.size(sim[[i]])-2)), L.base=~`+`, Ls.path=c(~`-`,~`-`)))/edges[i]
     sim_esf_neg[[i]] <- summary_formula(tmp ~ espL(c(0:(network.size(sim[[i]])-2)), L.base=~`-`, Ls.path=c(~`+`,~`+`)))/edges[i]
     sim_esf_pos[[i]] <-  summary_formula(tmp ~ espL(c(0:(network.size(sim[[i]])-2)), L.base=~`+`, Ls.path=c(~`+`,~`+`)))/edges[i]
   }
  sim_degree_pos <- bind_rows(sim_degree_pos)
  sim_degree_neg <- bind_rows(sim_degree_neg)
  sim_ese_neg <- bind_rows(sim_ese_neg)
  sim_ese_pos <- bind_rows(sim_ese_pos)
  sim_esf_neg <- bind_rows(sim_esf_neg)
  sim_esf_pos <- bind_rows(sim_esf_pos)

  colnames(sim_degree_pos) <- seq(0, ncol(sim_degree_pos))
  colnames(sim_degree_neg) <- seq(0, ncol(sim_degree_neg))
  colnames(sim_ese_neg) <- seq(0, ncol(sim_ese_neg))
  colnames(sim_ese_pos) <- seq(0, ncol(sim_ese_pos))
  colnames(sim_esf_neg) <- seq(0, ncol(sim_esf_neg))
  colnames(sim_esf_pos) <- seq(0, ncol(sim_esf_pos))

  # Extract max index
  get_max_indices <- function(x) {
    # Get column indices where the absolute sum is not zero
    col_indices <- which(sapply(x, function(y) sum(abs(y))) != 0)
    # Get the index of the maximum value
    max_indices <- max(col_indices)
    # Return the max indices
    return(max_indices)
  }


  # plot network statistics
  boxplot(sim_degree_pos,
          names = names(sim_degree_pos),
          xlab = "Positive degrees",
          ylab = "Proportion of nodes",
          ylim = c(0,max(unlist(degree_pos), unlist(sim_degree_pos), na.rm = T)),
          xlim = c(0, max(get_max_indices(sim_degree_pos), get_max_indices(degree_pos)) + 0.5))
  lines(degree_pos, lwd = 3)

  boxplot(sim_degree_neg,
          names = names(sim_degree_neg),
       xlab = "Negative degrees",
       ylab = "Proportion of nodes",
       ylim = c(0,max(unlist(sim_degree_neg), unlist(degree_neg), na.rm = T)),
       xlim = c(0, max(get_max_indices(sim_degree_neg), get_max_indices(degree_neg))+ 0.5))
  lines(degree_neg, lwd = 3)

  boxplot(sim_ese_neg,
          names = names(sim_ese_neg),
          xlab = "k−Edgewise-Shared Enemies (-)",
          ylab = "Proportion of edges",
          ylim = c(0, max(unlist(sim_ese_neg), unlist(ese_neg), na.rm = T)),
          xlim = c(0, max(get_max_indices(sim_ese_neg), get_max_indices(ese_neg))+ 0.5))
  lines(ese_neg, lwd = 3)

  boxplot(sim_ese_pos,
          names = names(sim_ese_pos),
          xlab = "k−Edgewise-Shared Enemies (+)",
          ylab = "Proportion of edges",
          ylim = c(0,max(unlist(sim_ese_pos), unlist(ese_pos), na.rm = T)),
          xlim = c(0, max(get_max_indices(sim_ese_pos), get_max_indices(ese_pos))+ 0.5))
  lines(ese_pos, lwd = 3)

  boxplot(sim_esf_neg,
          names = names(sim_esf_neg),
          xlab = "k−Edgewise-Shared Friends (-)",
          ylab = "Proportion of edges",
          ylim = c(0,max(unlist(sim_esf_neg), unlist(esf_neg), na.rm = T)),
          xlim = c(0, max(get_max_indices(sim_esf_neg), get_max_indices(esf_neg))+ 0.5))
  lines(esf_neg, lwd = 3)

  boxplot(sim_esf_pos,
          names = names(sim_esf_pos),
          xlab = "k−Edgewise-Shared Friends (+)",
          ylab = "Proportion of edges",
          ylim = c(0,max(unlist(sim_esf_pos), unlist(esf_pos), na.rm = T)),
          xlim = c(0, max(get_max_indices(sim_esf_pos), get_max_indices(esf_pos))+ 0.5))
  lines(esf_pos, lwd = 3)
}

