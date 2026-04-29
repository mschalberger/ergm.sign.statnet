#' Conduct Goodness-of-Fit Diagnostics for a Signed ERGM
#'
#' @param model A fitted signed ERGM (SERGM) object of class \code{"sign"}.
#' @param ... Additional arguments passed to methods.
#' @export
gof <- function(model, ...) {
  UseMethod("gof")
}

#' @rdname gof
#' @param nsim Integer; number of simulated networks. Defaults to 200.
#' @param seed Optional integer seed.
#' @return An object of class \code{"gof.sign"}.
#'   Print shows a summary; \code{plot()} produces diagnostic boxplots.
#' @export
gof.sign <- function(model, nsim = 200, seed = NULL, ...) {

  net      <- model$network
  dynamic  <- "dynamic.sign" %in% class(net)
  multi    <- "multi.sign"   %in% class(net)

  sim <- simulate(model,
                  nsim = nsim, seed = seed, basis = net)

  dispatch  <- .get_dispatch(net)
  #clean     <- dispatch$clean
  directed  <- .is_directed(net)

  obs      <- .stats_network(net, dynamic, multi, directed)
  sim_list <- lapply(seq_len(nsim), function(i)
    .stats_network(sim[[i]], dynamic, multi, directed)
  )

  stat_names <- names(obs)
  sim_df <- setNames(
    lapply(stat_names, function(nm)
      dplyr::bind_rows(lapply(sim_list, `[[`, nm))
    ),
    stat_names
  )

  structure(
    list(
      obs      = obs,
      sim      = sim_df,
      directed = directed,
      dynamic  = dynamic,
      multi    = multi,
      panels   = .build_panels(directed)
    ),
    class = "gof.sign"
  )
}

# =============================================================================
# S3 methods
# =============================================================================

#' @export
print.gof.sign <- function(x, ...) {
  cat("Goodness-of-fit diagnostics for a signed ERGM\n")
  cat("  Network type :", if (x$directed) "directed" else "undirected",
      if (x$dynamic) "(dynamic)" else if (x$multi) "(multi)" else "(static)", "\n")
  cat("  Statistics   :", paste(names(x$obs), collapse = ", "), "\n")
  cat("  Simulations  :", nrow(x$sim[[1]]), "\n")
  invisible(x)
}

#' @export
plot.gof.sign <- function(x, which = NULL, ...) {
  panels <- x$panels
  if (!is.null(which)) panels <- panels[which]
  for (p in panels)
    .plot_box(x$sim[[p$stat]], x$obs[[p$stat]], p$xlab, p$ylab)
  invisible(x)
}

# =============================================================================
# Network → stats: entry point handling static / dynamic / multi
# =============================================================================

#' @keywords internal
.stats_network <- function(net, dynamic, multi, directed) {
  if (dynamic || multi) {
    # UnLayer gives a list of static signed graphs (length = t-1 transitions)
    sgl    <- UnLayer(net)
    # remove first element
    sgl <- sgl[-1]
    A_list <- lapply(sgl, function(g) as.sociomatrix(g, attrname = "sign"))
    .stats_dynamic(A_list, directed)
  } else {
    A         <- as.sociomatrix(net, attrname = "sign")
    .stats_from_adj(A, directed)
  }
}

# =============================================================================
# Per-timepoint normalisation for dynamic / multi networks
# =============================================================================

#' @keywords internal
.stats_dynamic <- function(A_list, directed) {
  # Compute stats on each wave's adjacency matrix independently — each wave
  # uses its own edge count as denominator — then average the proportions.
  # Waves with few edges contribute equally to a wave with many, so density
  # differences across time do not distort the aggregate GOF statistic.
  wave_stats <- lapply(A_list, .stats_from_adj, directed = directed)

  stat_names <- names(wave_stats[[1]])
  setNames(
    lapply(stat_names, function(nm) {
      mat <- do.call(rbind, lapply(wave_stats, `[[`, nm))
      colSums(mat, na.rm = TRUE)
    }),
    stat_names
  )
}

# =============================================================================
# Matrix-based stat computation from a single signed adjacency matrix
# =============================================================================

#' @keywords internal
.stats_from_adj <- function(A, directed) {
  # A : n×n integer matrix with values in {-1, 0, +1}
  n    <- nrow(A)
  kmax <- n - 2L
  Ap   <- (A ==  1L) * 1L
  An   <- (A == -1L) * 1L
  esp  <- .compute_esp_stats(Ap, An, directed, kmax)

  if (directed) {
    list(
      idegree_pos = .degree_dist(Ap, n, "in"),
      idegree_neg = .degree_dist(An, n, "in"),
      odegree_pos = .degree_dist(Ap, n, "out"),
      odegree_neg = .degree_dist(An, n, "out"),
      esf_pos     = esp$esf_pos,
      esf_neg     = esp$esf_neg,
      ese_pos     = esp$ese_pos,
      ese_neg     = esp$ese_neg
    )
  } else {
    list(
      degree_pos  = .degree_dist(Ap, n, "total"),
      degree_neg  = .degree_dist(An, n, "total"),
      esf_pos     = esp$esf_pos,
      esf_neg     = esp$esf_neg,
      ese_pos     = esp$ese_pos,
      ese_neg     = esp$ese_neg
    )
  }
}

# =============================================================================
# Degree distributions
# =============================================================================

#' @keywords internal
.degree_dist <- function(Abin, n, type = c("total", "in", "out")) {
  type <- match.arg(type)
  deg  <- switch(type,
                 total = rowSums(Abin) + colSums(Abin),
                 out   = rowSums(Abin),
                 `in`  = colSums(Abin)
  )
  tab <- tabulate(deg + 1L, nbins = n)
  setNames(tab / n, 0:(n - 1L))
}

# =============================================================================
# ESP distributions
# =============================================================================

#' @keywords internal
.all_path_types <- function(Abin) {
  # Sum OTP + ITP + OSP + ISP, matching espL behaviour for directed networks
  At  <- t(Abin)
  OTP <- Abin %*% Abin   # i->k,  k->j
  ITP <- At   %*% At     # j->k,  k->i
  OSP <- Abin %*% At     # i->k,  j->k
  ISP <- At   %*% Abin   # k->i,  k->j
  OTP + ITP + OSP + ISP
}

#' @keywords internal
.esp_dist_mat <- function(closing_layer, path_mat, kmax) {
  ep <- sum(closing_layer)
  if (ep == 0L) return(setNames(rep(0, kmax + 1L), 0:kmax))
  counts <- pmin(path_mat[closing_layer == 1L], kmax)
  tab    <- tabulate(counts + 1L, nbins = kmax + 1L)
  setNames(tab / ep, 0:kmax)
}

#' @keywords internal
.compute_esp_stats <- function(Ap, An, directed, kmax) {
  PP <- if (directed) .all_path_types(Ap) else Ap %*% Ap
  NN <- if (directed) .all_path_types(An) else An %*% An
  list(
    esf_pos = .esp_dist_mat(Ap, PP, kmax),   # shared friends of + edge
    esf_neg = .esp_dist_mat(An, PP, kmax),   # shared friends of - edge
    ese_pos = .esp_dist_mat(Ap, NN, kmax),   # shared enemies of + edge
    ese_neg = .esp_dist_mat(An, NN, kmax)    # shared enemies of - edge
  )
}

# =============================================================================
# Dispatch
# =============================================================================

#' @keywords internal
.get_dispatch <- function(net) {
  if ("static.sign" %in% class(net)) {
    list(clean = identity)
  } else if ("dynamic.sign" %in% class(net)) {
    list(clean = function(x) { class(x) <- setdiff(class(x), "dynamic.sign"); x })
  } else if ("multi.sign" %in% class(net)) {
    list(clean = function(x) { class(x) <- setdiff(class(x), "multi.sign"); x })
  } else {
    stop("Unsupported network class.")
  }
}

#' @keywords internal
.is_directed <- function(net) isTRUE(net %n% "directed")

# =============================================================================
# Panel metadata
# =============================================================================

#' @keywords internal
.build_panels <- function(directed) {
  if (directed) {
    list(
      list(stat = "idegree_pos", xlab = "Positive in-degree",        ylab = "Proportion of nodes"),
      list(stat = "idegree_neg", xlab = "Negative in-degree",        ylab = "Proportion of nodes"),
      list(stat = "odegree_pos", xlab = "Positive out-degree",       ylab = "Proportion of nodes"),
      list(stat = "odegree_neg", xlab = "Negative out-degree",       ylab = "Proportion of nodes"),
      list(stat = "esf_pos",     xlab = "k-ESF: shared friends (+)", ylab = "Proportion of positive edges"),
      list(stat = "esf_neg",     xlab = "k-ESF: shared friends (-)", ylab = "Proportion of negative edges"),
      list(stat = "ese_pos",     xlab = "k-ESE: shared enemies (+)", ylab = "Proportion of positive edges"),
      list(stat = "ese_neg",     xlab = "k-ESE: shared enemies (-)", ylab = "Proportion of negative edges")
    )
  } else {
    list(
      list(stat = "degree_pos",  xlab = "Positive degree",           ylab = "Proportion of nodes"),
      list(stat = "degree_neg",  xlab = "Negative degree",           ylab = "Proportion of nodes"),
      list(stat = "esf_pos",     xlab = "k-ESF: shared friends (+)", ylab = "Proportion of positive edges"),
      list(stat = "esf_neg",     xlab = "k-ESF: shared friends (-)", ylab = "Proportion of negative edges"),
      list(stat = "ese_pos",     xlab = "k-ESE: shared enemies (+)", ylab = "Proportion of positive edges"),
      list(stat = "ese_neg",     xlab = "k-ESE: shared enemies (-)", ylab = "Proportion of negative edges")
    )
  }
}

# =============================================================================
# Plotting
# =============================================================================

#' @keywords internal
.plot_box <- function(sim_data, obs_data, xlab, ylab) {

  get_max_nonzero <- function(x) {
    cols <- which(vapply(as.data.frame(x),
                         function(y) sum(abs(y), na.rm = TRUE),
                         numeric(1)) != 0)
    if (length(cols) == 0L) 1L else max(cols)
  }

  ylim_max <- max(unlist(sim_data), unlist(obs_data), na.rm = TRUE)
  xlim_max <- max(get_max_nonzero(sim_data),
                  get_max_nonzero(obs_data)) + 0.5

  par(bty = "l")
  boxplot(sim_data,
          names = names(sim_data),
          xlab  = xlab,
          ylab  = ylab,
          ylim  = c(0, ylim_max),
          xlim  = c(0, xlim_max),
          xaxt  = "n",
          yaxt  = "n")
  axis(1, at = seq_along(sim_data), labels = 0:(length(sim_data) - 1L))
  axis(2, at = pretty(c(0, ylim_max)), las = 1)
  lines(obs_data, lwd = 3)
}
