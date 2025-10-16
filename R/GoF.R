#' Conduct Goodness-of-Fit Diagnostics on a Signed Exponential Family Random Graph Model
#'
#' This function calculates the goodness-of-fit (GoF) for a fitted SERGM model by simulating new networks using the same model and comparing the statistics of the original network with those of the simulated networks. The statistics used in the comparison are the positive and negative degree distributions, edge-wise shared enemies distributions for positive and negative edges, and edge-wise shared friend distributions for positive and negative edges.
#'
#' @param model A fitted SERGM model.
#' @param nsim An integer representing the number of simulated networks to generate. Defaults to 200.
#' @param seed Seed value (interger) for the random number generator. See \link{set.seed}
#'
#' @return Plots 6 diagnostics for the goodness-of-fit of signed exponential family random graph models.
#'
#' @seealso \link{sergm}, \link{gof_tsergm}
#'
#' @examples
#'data("tribes")
#'fit <- ergm(tribes ~ Pos(~ edges) + Neg(~ edges))
#'GoF(fit, nsim = 20)
#'
#' @export

GoF <- function(model, nsim = 200, seed = NULL) {

  # simulate
  sim <- simulate(model, nsim = nsim, seed = seed)

  if("static.sign" %in% class(sim[[1]])) {
    # get underlying network
    pos <- as.sociomatrix(model[["network"]][["gal"]][[".subnetcache"]][[".LayerID"]][["+"]])
    neg <- as.sociomatrix(model[["network"]][["gal"]][[".subnetcache"]][[".LayerID"]][["-"]])*-1
    comb <- pos + neg
    net <- network.sign(mat = comb, matrix.type = "adjacency")

    n <- ncol(pos)
    e <- summary_formula(net ~ edges)

    # create network statistics
    degree_pos <- summary_formula(net ~ Pos(~ degree(0:(n-1)))) / n
    degree_neg <- summary_formula(net ~ Neg(~ degree(0:(n-1)))) / n
    names(degree_neg) <- names(degree_pos) <- 0:(n-1)

    ese_neg <- summary_formula(net ~ ese(0:(n-2), base = "-")) / e
    ese_pos <- summary_formula(net ~ ese(0:(n-2), base="+")) / e
    esf_neg <- summary_formula(net ~ esf(0:(n-2), base="-")) / e
    esf_pos <- summary_formula(net ~ esf(0:(n-2), base="+")) / e
    names(esf_neg) <- names(ese_pos) <- names(ese_neg) <- names(esf_pos) <- 0:(n-2)

    # create network statistics for simulated networks
    sim_degree_pos <- sim_degree_neg <- net_size <- edges <- sim_ese_neg <- sim_ese_pos <- sim_esf_neg <- sim_esf_pos <- c()
    for (i in 1:nsim) {
      tmp <- sim[[i]]
      net_size[i] <- network.size(tmp)/2
      edges[i] <- summary_formula(tmp ~ edges)
      sim_degree_pos[[i]] <- summary_formula(tmp ~ Pos(~ degree(0:(net_size[i]-1))))/net_size[i]
      sim_degree_neg[[i]] <- summary_formula(tmp ~ Neg(~ degree(0:(net_size[i]-1))))/net_size[i]
      sim_ese_neg[[i]] <- summary_formula(tmp ~ ese(0:(net_size[i]-2), base="-"))/edges[i]
      sim_ese_pos[[i]] <- summary_formula(tmp ~ ese(0:(net_size[i]-2), base="+"))/edges[i]
      sim_esf_neg[[i]] <- summary_formula(tmp ~ esf(0:(net_size[i]-2), base="-"))/edges[i]
      sim_esf_pos[[i]] <-  summary_formula(tmp ~ esf(0:(net_size[i]-2), base="+"))/edges[i]
    }
    sim_degree_pos <- bind_rows(sim_degree_pos)
    sim_degree_neg <- bind_rows(sim_degree_neg)
    sim_ese_neg <- bind_rows(sim_ese_neg)
    sim_ese_pos <- bind_rows(sim_ese_pos)
    sim_esf_neg <- bind_rows(sim_esf_neg)
    sim_esf_pos <- bind_rows(sim_esf_pos)

    colnames(sim_degree_pos) <- colnames(sim_degree_neg) <- seq(0, ncol(sim_degree_neg))
    colnames(sim_ese_neg) <- colnames(sim_ese_pos) <- colnames(sim_esf_neg) <- colnames(sim_esf_pos) <- seq(0, ncol(sim_esf_pos))

    # Extract max index
    get_max_indices <- function(x) {
      col_indices <- which(sapply(x, function(y) sum(abs(y))) != 0)
      max_indices <- max(col_indices)
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
   else if ("dynamic.sign" %in% class(sim[[1]])) {
     net_list <- model[["network"]][["gal"]][["NetList"]]

    n <- e <- degree_pos <- degree_neg <- ese_neg <- ese_pos <- esf_neg <- esf_pos <- c()
    for (i in 1:length(net_list)) {
      MultiNet <- net_list[[i]]
      nw <- UnLayer(MultiNet)
      n[[i]] <- network.size(nw)

      e[[i]] <- summary_formula(MultiNet ~ edges)

      degree_pos[[i]] <- summary_formula(MultiNet ~ Pos(~ degree(0:(n[[i]]-1))))
      degree_neg[[i]] <- summary_formula(MultiNet ~ Neg(~ degree(0:(n[[i]]-1))))
      names(degree_pos[[i]]) <- names(degree_neg[[i]]) <- 0:(n[[i]]-1)

      ese_neg[[i]] <- summary_formula(MultiNet ~ ese(0:(n[[i]]-2), base="-"))
      ese_pos[[i]] <- summary_formula(MultiNet ~ ese(0:(n[[i]]-2), base="+"))
      esf_neg[[i]] <- summary_formula(MultiNet ~ esf(0:(n[[i]]-2), base="-"))
      esf_pos[[i]] <- summary_formula(MultiNet ~ esf(0:(n[[i]]-2), base="+"))
      names(esf_neg[[i]]) <- names(ese_pos[[i]]) <- names(ese_neg[[i]]) <- names(esf_pos[[i]]) <- 0:(n[[i]]-2)
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
      sim_list <- UnLayer(sim[[i]])
      sim_degree_pos[[i]] <- sim_degree_neg[[i]] <- net_size[[i]] <- edges[[i]] <- sim_ese_neg[[i]] <- sim_ese_pos[[i]] <- sim_esf_neg[[i]] <- sim_esf_pos[[i]] <- list()
      for (j in 1:length(sim_list$multi)) {
        net_size[[i]][[j]] <- network.size(sim_list$single[[j]])
        tmp <- sim_list$multi[[j]]
        edges[[i]][[j]] <- summary_formula(tmp ~ edges)
        sim_degree_pos[[i]][[j]] <- summary_formula(tmp ~ Pos(~ degree(0:(net_size[[i]][[j]]-1))))
        sim_degree_neg[[i]][[j]] <- summary_formula(tmp ~ Neg(~ degree(0:(net_size[[i]][[j]]-1))))
        sim_ese_neg[[i]][[j]] <- summary_formula(tmp ~ ese(0:(net_size[[i]][[j]]-2), base="-"))
        sim_ese_pos[[i]][[j]] <- summary_formula(tmp ~ ese(0:(net_size[[i]][[j]]-2), base="+"))
        sim_esf_neg[[i]][[j]] <- summary_formula(tmp ~ esf(0:(net_size[[i]][[j]]-2), base="-"))
        sim_esf_pos[[i]][[j]] <-  summary_formula(tmp ~ esf(0:(net_size[[i]][[j]]-2), base="+"))
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
      col_indices <- which(sapply(x, function(y) sum(abs(y))) != 0)
      max_indices <- max(col_indices)
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
   }
}


