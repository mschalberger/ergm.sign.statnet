#' Pool Timepoints of a Dynamic Signed Network
#'
#' Aggregates timepoints of a dynamic signed network into pooled periods using
#' majority-vote logic, with configurable tie-breaking.
#'
#' @param net A dynamic signed network of class `dynamic.sign`.
#' @param timepoints Either an integer specifying the number of equal-sized groups,
#'   or a list of integer vectors specifying which timepoints belong to each group.
#'   E.g., `list(c(1,2,3), c(4,5,6))` or `3` for three equal groups.
#' @param tie.breaker How to resolve ties when positive and negative counts are equal.
#'   One of `"zero"` (default), `"positive"`, `"negative"`, `"first"`, `"last"`.
#' @param dual.sign Logical. Allow positive and negative edges simultaneously
#'   between the same pair. Defaults to FALSE.
#'
#' @return A dynamic signed network of class `dynamic.sign` with one network per
#'   pooled group, or a `static.sign` network if only one group results.
#' @export
pool.sign <- function(net,
                      timepoints,
                      tie.breaker = c("zero", "positive", "negative", "first", "last"),
                      dual.sign = FALSE) {

  tie.breaker <- match.arg(tie.breaker)

  # --- Input validation ---
  if (!inherits(net, "dynamic.sign")) {
    stop("`net` must be a dynamic signed network of class 'dynamic.sign'.")
  }

  nets         <- UnLayer(net)
  n_timepoints <- length(nets)

  if (n_timepoints < 2) {
    stop("`net` must contain at least two timepoints to pool.")
  }

  # --- Extract signed adjacency matrices (1 = pos, -1 = neg, 2 = both) ---
  mats <- lapply(seq_len(n_timepoints), function(i) {
    as.sociomatrix(nets[[i]], attrname = "sign")
  })

  to_pos <- function(m) (m ==  1 | m == 2) * 1
  to_neg <- function(m) (m == -1 | m == 2) * 1

  pos_mats <- lapply(mats, to_pos)
  neg_mats <- lapply(mats, to_neg)

  # --- Resolve timepoints argument ---
  if (is.numeric(timepoints) && length(timepoints) == 1) {
    n_groups   <- as.integer(timepoints)
    group_size <- floor(n_timepoints / n_groups)
    timepoints <- split(seq_len(n_timepoints),
                        ceiling(seq_len(n_timepoints) / group_size))
  } else if (is.list(timepoints)) {
    all_idx <- unlist(timepoints)
    if (any(all_idx < 1) || any(all_idx > n_timepoints)) {
      stop("All timepoint indices must be between 1 and ", n_timepoints, ".")
    }
  } else {
    stop("`timepoints` must be a single integer (number of groups) or a list of index vectors.")
  }

  # --- Majority-vote pooling over a group of matrices ---
  pool_group <- function(mats, group_indices, tie.breaker, dual.sign = FALSE) {
    group_mats <- mats[group_indices]
    mat_dim    <- dim(group_mats[[1]])
    mat_array  <- array(unlist(group_mats), dim = c(mat_dim[1], mat_dim[2], length(group_mats)))

    pos_count <- apply(mat_array, 1:2, function(x) sum(x > 0))
    neg_count <- apply(mat_array, 1:2, function(x) sum(x < 0))

    result <- matrix(0, nrow = mat_dim[1], ncol = mat_dim[2])
    result[pos_count > neg_count] <-  1
    result[neg_count > pos_count] <- -1

    tie_mask <- (pos_count == neg_count) & (pos_count > 0)
    if (any(tie_mask)) {
      if (dual.sign) {
        # When dual.sign is active, ties become dual edges (both positive and negative)
        result[tie_mask] <- 2L
      } else {
        result[tie_mask] <- switch(tie.breaker,
                                   zero     = 0,
                                   positive =  1,
                                   negative = -1,
                                   first    = group_mats[[1]][tie_mask],
                                   last     = group_mats[[length(group_mats)]][tie_mask],
                                   stop("Unknown tie.breaker: ", tie.breaker)
        )
      }
    }

    ref_dn <- dimnames(group_mats[[1]])
    if (!is.null(ref_dn[[1]]) && !is.null(ref_dn[[2]])) {
      dimnames(result) <- ref_dn
    }

    attr(result, "pool_label") <- paste0("Timepoints ", min(group_indices), "-", max(group_indices))
    result
  }

  pooled_pos <- lapply(timepoints, function(g) pool_group(pos_mats, g, tie.breaker))
  pooled_neg <- lapply(timepoints, function(g) pool_group(neg_mats, g, tie.breaker))

  if (length(pooled_pos) == 1) {
    # If only one group, return a static.sign network
    result_net <- network.sign(
      pos.mat   = pooled_pos[[1]],
      neg.mat   = pooled_neg[[1]],
      directed  = net %n% "directed",
      loops     = net %n% "loops",
      dual.sign = dual.sign
    )

    result_net %n% "names"     <- attr(pooled_pos[[1]], "pool_label")
    result_net %n% "dual.sign" <- dual.sign
  } else {
  result_net <- network.sign(
    pos.mat   = pooled_pos,
    neg.mat   = pooled_neg,
    directed  = net %n% "directed",
    loops     = net %n% "loops",
    dual.sign = dual.sign
  )

  result_net %n% "names"     <- sapply(pooled_pos, attr, "pool_label")
  result_net %n% "dual.sign" <- dual.sign
  }

  return(result_net)
}
