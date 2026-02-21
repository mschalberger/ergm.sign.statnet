#' Conduct Goodness-of-Fit Diagnostics for a Signed ERGM
#'
#' Computes the goodness-of-fit (GOF) for a fitted signed exponential random
#' graph model (SERGM). The function simulates new networks using the fitted
#' model and compares key network statistics from the observed network with
#' those from the simulated ones.
#'
#' The following diagnostics are plotted:
#' \itemize{
#'   \item Positive degree distribution
#'   \item Negative degree distribution
#'   \item Edgewise shared enemies distribution (positive edges)
#'   \item Edgewise shared enemies distribution (negative edges)
#'   \item Edgewise shared friends distribution (positive edges)
#'   \item Edgewise shared friends distribution (negative edges)
#' }
#'
#' @param model A fitted signed ERGM (SERGM) object.
#' @param nsim Integer; number of simulated networks to generate. Defaults to 200.
#' @param seed Optional integer seed for reproducibility. Passed to
#'   \code{\link{set.seed}}.
#'
#' @return Produces six diagnostic boxplots comparing observed and simulated
#'   statistics for the fitted model.
#'
#' @seealso \code{\link[ergm]{ergm}}, \code{\link{mple_sign}}
#'
#' @examples
#'\donttest{
#' data("tribes")
#' fit <- mple_sign(tribes ~ Pos(~edges) + Neg(~edges))
#' gof.sign(fit, nsim = 100)
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom stats simulate
#' @importFrom graphics axis
#' @export
gof.sign <- function(model, nsim = 200, seed = NULL) {

  # simulate networks
  net <- model$network
  sim <- simulate(model$formula, coef = unname(model$coefficients), nsim = nsim, seed = seed, basis = net,
                  constraints = net[["gal"]][["ergm"]][["constraints"]])

  # select appropriate formula based on model type
  if ("static.sign" %in% class(net)) {
    n <- network.size(net) / 2
    #e <- summary_formula(net ~ edges)
    e_pos <- summary_formula(net ~ Pos(~edges))
    e_neg <- summary_formula(net ~ Neg(~edges))

    degree_pos <- summary_formula(net ~ Pos(~degree(0:(n - 1)))) / n
    degree_neg <- summary_formula(net ~ Neg(~degree(0:(n - 1)))) / n
    names(degree_neg) <- names(degree_pos) <- 0:(n - 1)

    ese_neg <- summary_formula(net ~ ese(0:(n - 2), base = "-")) / e_neg
    ese_pos <- summary_formula(net ~ ese(0:(n - 2), base = "+")) / e_pos
    esf_neg <- summary_formula(net ~ esf(0:(n - 2), base = "-")) / e_neg
    esf_pos <- summary_formula(net ~ esf(0:(n - 2), base = "+")) / e_pos
    names(esf_neg) <- names(ese_pos) <- names(ese_neg) <- names(esf_pos) <- 0:(n - 2)

    sim_degree_pos <- sim_degree_neg <- sim_ese_neg <- sim_ese_pos <- sim_esf_neg <- sim_esf_pos <- vector("list", nsim)
    edges <- edges_pos <- edges_neg <- numeric(nsim)

    for (i in seq_len(nsim)) {
      tmp <- sim[[i]]
      #edges[i] <- summary_formula(tmp ~ edges)
      edges_pos[i] <- summary_formula(tmp ~ Pos(~edges))
      edges_neg[i] <- summary_formula(tmp ~ Neg(~edges))
      sim_degree_pos[[i]] <- summary_formula(tmp ~ Pos(~degree(0:(n - 1)))) / n
      sim_degree_neg[[i]] <- summary_formula(tmp ~ Neg(~degree(0:(n - 1)))) / n
      sim_ese_neg[[i]] <- summary_formula(tmp ~ ese(0:(n - 2), base = "-")) / edges_neg[i]
      sim_ese_pos[[i]] <- summary_formula(tmp ~ ese(0:(n - 2), base = "+")) / edges_pos[i]
      sim_esf_neg[[i]] <- summary_formula(tmp ~ esf(0:(n - 2), base = "-")) / edges_neg[i]
      sim_esf_pos[[i]] <- summary_formula(tmp ~ esf(0:(n - 2), base = "+")) / edges_pos[i]
    }

  } else if ("dynamic.sign" %in% class(net)) {
    class(net) <- setdiff(class(net), "dynamic.sign")

    n <- network.size(net) / 2
    #e <- summary_formula(net ~ Cross(~edges))
    e_pos <- summary_formula(net ~ Cross(~Pos(~edges)))
    e_neg <- summary_formula(net ~ Cross(~Neg(~edges)))

    degree_pos <- summary_formula(net ~ Cross(~Pos(~degree(0:(n - 1))))) / n
    degree_neg <- summary_formula(net ~ Cross(~Neg(~degree(0:(n - 1))))) / n
    names(degree_neg) <- names(degree_pos) <- 0:(n - 1)

    ese_neg <- summary_formula(net ~ Cross(~ese(0:(n - 2), base = "-"))) / e_neg
    ese_pos <- summary_formula(net ~ Cross(~ese(0:(n - 2), base = "+"))) / e_pos
    esf_neg <- summary_formula(net ~ Cross(~esf(0:(n - 2), base = "-"))) / e_neg
    esf_pos <- summary_formula(net ~ Cross(~esf(0:(n - 2), base = "+"))) / e_pos
    names(esf_neg) <- names(ese_pos) <- names(ese_neg) <- names(esf_pos) <- 0:(n - 2)

    sim_degree_pos <- sim_degree_neg <- sim_ese_neg <- sim_ese_pos <- sim_esf_neg <- sim_esf_pos <- vector("list", nsim)
    edges <- edges_pos <- edges_neg <- numeric(nsim)

    for (i in seq_len(nsim)) {
      tmp <- sim[[i]]
      class(tmp) <- setdiff(class(tmp), "dynamic.sign")
      #edges[i] <- summary_formula(tmp ~ Cross(~edges))
      edges_pos[i] <- summary_formula(tmp ~ Cross(~Pos(~edges)))
      edges_neg[i] <- summary_formula(tmp ~ Cross(~Neg(~edges)))
      sim_degree_pos[[i]] <- summary_formula(tmp ~ Cross(~Pos(~degree(0:(n - 1))))) / n
      sim_degree_neg[[i]] <- summary_formula(tmp ~ Cross(~Neg(~degree(0:(n - 1))))) / n
      sim_ese_neg[[i]] <- summary_formula(tmp ~ Cross(~ese(0:(n - 2), base = "-"))) / edges_neg[i]
      sim_ese_pos[[i]] <- summary_formula(tmp ~ Cross(~ese(0:(n - 2), base = "+"))) / edges_pos[i]
      sim_esf_neg[[i]] <- summary_formula(tmp ~ Cross(~esf(0:(n - 2), base = "-"))) / edges_neg[i]
      sim_esf_pos[[i]] <- summary_formula(tmp ~ Cross(~esf(0:(n - 2), base = "+"))) / edges_pos[i]
    }

  } else if ("multi.sign" %in% class(net)) {
    class(net) <- setdiff(class(net), "multi.sign")

    n <- network.size(net) / 2
    #e <- summary_formula(net ~ N(~edges))
    e_pos <- summary_formula(net ~ N(~Pos(~edges)))
    e_neg <- summary_formula(net ~ N(~Neg(~edges)))

    degree_pos <- summary_formula(net ~ N(~Pos(~degree(0:(n - 1))))) / n
    degree_neg <- summary_formula(net ~ N(~Neg(~degree(0:(n - 1))))) / n
    names(degree_neg) <- names(degree_pos) <- 0:(n - 1)

    ese_neg <- summary_formula(net ~ N(~ese(0:(n - 2), base = "-"))) / e_neg
    ese_pos <- summary_formula(net ~ N(~ese(0:(n - 2), base = "+"))) / e_pos
    esf_neg <- summary_formula(net ~ N(~esf(0:(n - 2), base = "-"))) / e_neg
    esf_pos <- summary_formula(net ~ N(~esf(0:(n - 2), base = "+"))) / e_pos
    names(esf_neg) <- names(ese_pos) <- names(ese_neg) <- names(esf_pos) <- 0:(n - 2)

    sim_degree_pos <- sim_degree_neg <- sim_ese_neg <- sim_ese_pos <- sim_esf_neg <- sim_esf_pos <- vector("list", nsim)
    edges <- edges_pos <- edges_neg <- numeric(nsim)

    for (i in seq_len(nsim)) {
      tmp <- sim[[i]]
      class(tmp) <- setdiff(class(tmp), "multi.sign")
      #edges[i] <- summary_formula(tmp ~ N(~edges))
      edges_pos[i] <- summary_formula(tmp ~ N(~Pos(~edges)))
      edges_neg[i] <- summary_formula(tmp ~ N(~Neg(~edges)))
      sim_degree_pos[[i]] <- summary_formula(tmp ~ N(~Pos(~degree(0:(n - 1))))) / n
      sim_degree_neg[[i]] <- summary_formula(tmp ~ N(~Neg(~degree(0:(n - 1))))) / n
      sim_ese_neg[[i]] <- summary_formula(tmp ~ N(~ese(0:(n - 2), base = "-"))) / edges_neg[i]
      sim_ese_pos[[i]] <- summary_formula(tmp ~ N(~ese(0:(n - 2), base = "+"))) / edges_pos[i]
      sim_esf_neg[[i]] <- summary_formula(tmp ~ N(~esf(0:(n - 2), base = "-"))) / edges_neg[i]
      sim_esf_pos[[i]] <- summary_formula(tmp ~ N(~esf(0:(n - 2), base = "+"))) / edges_pos[i]
    }
  } else {
    stop("Unsupported network class.")
  }

  # combine simulation results
  sim_degree_pos <- bind_rows(sim_degree_pos)
  sim_degree_neg <- bind_rows(sim_degree_neg)
  sim_ese_neg <- bind_rows(sim_ese_neg)
  sim_ese_pos <- bind_rows(sim_ese_pos)
  sim_esf_neg <- bind_rows(sim_esf_neg)
  sim_esf_pos <- bind_rows(sim_esf_pos)

  # helper to determine nonzero columns
  get_max_indices <- function(x) {
    cols <- which(sapply(x, function(y) sum(abs(y))) != 0)
    max(cols)
  }

  # plotting
  plot_box <- function(sim_data, obs_data, xlab, ylab) {
    boxplot(sim_data,
            names = names(sim_data),
            xlab = xlab,
            ylab = ylab,
            ylim = c(0, max(unlist(sim_data), unlist(obs_data), na.rm = TRUE)),
            xlim = c(0, max(get_max_indices(sim_data), get_max_indices(obs_data)) + 0.5),
            xaxt = "n")               # suppress x-axis
    # add numeric x-axis from 0 to n
    axis(1, at = seq_along(sim_data), labels = 0:(length(sim_data) - 1))
    lines(obs_data, lwd = 3)
  }

  plot_box(sim_degree_pos, degree_pos, "Positive degrees", "Proportion of nodes")
  plot_box(sim_degree_neg, degree_neg, "Negative degrees", "Proportion of nodes")
  plot_box(sim_ese_neg, ese_neg, "k-Edgewise Shared Enemies (-)", "Proportion of negative edges")
  plot_box(sim_ese_pos, ese_pos, "k-Edgewise Shared Enemies (+)", "Proportion of positive edges")
  plot_box(sim_esf_neg, esf_neg, "k-Edgewise Shared Friends (-)", "Proportion of negative edges")
  plot_box(sim_esf_pos, esf_pos, "k-Edgewise Shared Friends (+)", "Proportion of positive edges")
}



