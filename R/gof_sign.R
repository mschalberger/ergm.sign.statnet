#' Conduct Goodness-of-Fit Diagnostics for a Signed ERGM
#'
#' @param model A fitted signed ERGM (SERGM) object of class \code{"sign"}.
#' @param ... Additional arguments passed to methods.
#' @export
gof <- function(model, ...) {
  UseMethod("gof")
}

#' @rdname gof
#' @param nsim Integer; number of simulated networks to generate. Defaults to 200.
#' @param seed Optional integer seed for reproducibility. Passed to \code{\link{set.seed}}.
#' @return Produces six diagnostic boxplots comparing observed and simulated
#'   statistics for the fitted model.
#' @seealso \code{\link[ergm]{ergm}}, \code{\link{mple_sign}}
#' @examples
#' \donttest{
#' data("tribes")
#' fit <- mple_sign(tribes ~ Pos(~edges) + Neg(~edges))
#' gof(fit, nsim = 100)
#' }
#' @importFrom dplyr bind_rows
#' @importFrom stats simulate
#' @importFrom graphics axis
#' @export
gof.sign <- function(model, nsim = 200, seed = NULL, ...) {

  net <- model$network
  sim <- simulate(model$formula, coef = unname(model$coefficients),
                  nsim = nsim, seed = seed, basis = net)

  stats <- .compute_gof_stats(net, sim, nsim)

  .plot_gof(stats)
}

# --- Internal helpers ---------------------------------------------------------

#' @keywords internal
.get_dispatch <- function(net) {
  if ("static.sign" %in% class(net)) {
    list(
      clean = identity,
      wrap  = function(term) term
    )
  } else if ("dynamic.sign" %in% class(net)) {
    list(
      clean = function(x) { class(x) <- setdiff(class(x), "dynamic.sign"); x },
      wrap  = function(term) bquote(Cross(~.(term)))
    )
  } else if ("multi.sign" %in% class(net)) {
    list(
      clean = function(x) { class(x) <- setdiff(class(x), "multi.sign"); x },
      wrap  = function(term) bquote(N(~.(term)))
    )
  } else {
    stop("Unsupported network class.")
  }
}

#' @keywords internal
.sf <- function(net, rhs) {
  summary_formula(as.formula(bquote(net ~ .(rhs)), env = environment()))
}

#' @keywords internal
.observed_stats <- function(net, n, wrap) {
  ep <- .sf(net, wrap(quote(Pos(~edges))))
  en <- .sf(net, wrap(quote(Neg(~edges))))

  deg_pos <- .sf(net, wrap(quote(Pos(~degree(0:(n - 1))))))
  deg_neg <- .sf(net, wrap(quote(Neg(~degree(0:(n - 1))))))

  # derive n from the actual output length rather than network.size
  n_deg <- length(deg_pos)
  deg_pos <- deg_pos / n
  deg_neg <- deg_neg / n
  names(deg_pos) <- names(deg_neg) <- 0:(n_deg - 1)

  ese_pos <- .sf(net, wrap(quote(ese(0:(n - 2), base = "+")))) / ep
  ese_neg <- .sf(net, wrap(quote(ese(0:(n - 2), base = "-")))) / en
  esf_pos <- .sf(net, wrap(quote(esf(0:(n - 2), base = "+")))) / ep
  esf_neg <- .sf(net, wrap(quote(esf(0:(n - 2), base = "-")))) / en

  n_ese <- length(ese_pos)
  names(ese_pos) <- names(ese_neg) <- names(esf_pos) <- names(esf_neg) <- 0:(n_ese - 1)

  list(degree_pos = deg_pos, degree_neg = deg_neg,
       ese_pos = ese_pos,    ese_neg = ese_neg,
       esf_pos = esf_pos,    esf_neg = esf_neg)
}

#' @keywords internal
.compute_gof_stats <- function(net, sim, nsim) {

  dispatch  <- .get_dispatch(net)
  clean     <- dispatch$clean
  wrap      <- dispatch$wrap

  net_clean <- clean(net)
  n         <- network.size(net_clean) / 2

  obs       <- .observed_stats(net_clean, n, wrap)
  sim_lists <- .empty_sim_lists(nsim)

  for (i in seq_len(nsim)) {
    tmp <- clean(sim[[i]])
    ep  <- .sf(tmp, wrap(quote(Pos(~edges))))
    en  <- .sf(tmp, wrap(quote(Neg(~edges))))

    dp <- .sf(tmp, wrap(quote(Pos(~degree(0:(n - 1))))))
    dn <- .sf(tmp, wrap(quote(Neg(~degree(0:(n - 1))))))
    sim_lists$degree_pos[[i]] <- setNames(dp / n, 0:(length(dp) - 1))
    sim_lists$degree_neg[[i]] <- setNames(dn / n, 0:(length(dn) - 1))

    ep_ese <- .sf(tmp, wrap(quote(ese(0:(n - 2), base = "+"))))
    en_ese <- .sf(tmp, wrap(quote(ese(0:(n - 2), base = "-"))))
    ep_esf <- .sf(tmp, wrap(quote(esf(0:(n - 2), base = "+"))))
    en_esf <- .sf(tmp, wrap(quote(esf(0:(n - 2), base = "-"))))
    sim_lists$ese_pos[[i]] <- setNames(ep_ese / ep, 0:(length(ep_ese) - 1))
    sim_lists$ese_neg[[i]] <- setNames(en_ese / en, 0:(length(en_ese) - 1))
    sim_lists$esf_pos[[i]] <- setNames(ep_esf / ep, 0:(length(ep_esf) - 1))
    sim_lists$esf_neg[[i]] <- setNames(en_esf / en, 0:(length(en_esf) - 1))
  }

  list(
    obs = obs,
    sim = lapply(sim_lists, bind_rows)
  )
}
#' @keywords internal
.empty_sim_lists <- function(nsim) {
  nms <- c("degree_pos", "degree_neg", "ese_neg", "ese_pos", "esf_neg", "esf_pos")
  setNames(replicate(length(nms), vector("list", nsim), simplify = FALSE), nms)
}

#' @keywords internal
.plot_gof <- function(stats) {
  panels <- list(
    list(sim = stats$sim$degree_pos, obs = stats$obs$degree_pos,
         xlab = "Positive degrees",              ylab = "Proportion of nodes"),
    list(sim = stats$sim$degree_neg, obs = stats$obs$degree_neg,
         xlab = "Negative degrees",              ylab = "Proportion of nodes"),
    list(sim = stats$sim$ese_neg,    obs = stats$obs$ese_neg,
         xlab = "k-Edgewise Shared Enemies (-)", ylab = "Proportion of negative edges"),
    list(sim = stats$sim$ese_pos,    obs = stats$obs$ese_pos,
         xlab = "k-Edgewise Shared Enemies (+)", ylab = "Proportion of positive edges"),
    list(sim = stats$sim$esf_neg,    obs = stats$obs$esf_neg,
         xlab = "k-Edgewise Shared Friends (-)", ylab = "Proportion of negative edges"),
    list(sim = stats$sim$esf_pos,    obs = stats$obs$esf_pos,
         xlab = "k-Edgewise Shared Friends (+)", ylab = "Proportion of positive edges")
  )

  for (p in panels) .plot_box(p$sim, p$obs, p$xlab, p$ylab)
}

#' @keywords internal
.plot_box <- function(sim_data, obs_data, xlab, ylab) {
  get_max_nonzero <- function(x) {
    cols <- which(vapply(x, function(y) sum(abs(y)), numeric(1)) != 0)
    max(cols)
  }

  par(bty = "l")
  boxplot(sim_data,
          names  = names(sim_data),
          xlab   = xlab,
          ylab   = ylab,
          ylim   = c(0, max(unlist(sim_data), unlist(obs_data), na.rm = TRUE)),
          xlim   = c(0, max(get_max_nonzero(sim_data),
                            get_max_nonzero(obs_data)) + 0.5),
          xaxt   = "n",
          yaxt   = "n")
  axis(1, at = seq_along(sim_data), labels = 0:(length(sim_data) - 1))
  axis(2, at = pretty(c(0, max(unlist(sim_data), unlist(obs_data), na.rm = TRUE))),
       las = 1)
  lines(obs_data, lwd = 3)
}



