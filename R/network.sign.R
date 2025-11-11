#' Create Signed Network Object
#'
#' Turn adjacency matrices or edgelists into static or dynamic signed networks.
#'
#' @param mat (List of) signed adjacency matrices or edgelists. For dynamic networks, provide a list.
#' Adjacency matrices must contain only -1, 0, or 1. Edgelists must have three columns: "From", "To", and "Sign".
#' @param pos.mat Optional. Positive adjacency matrix or list of matrices.
#' @param neg.mat Optional. Negative adjacency matrix or list of matrices.
#' If provided, these are treated as two separate layers of the same network.
#' @param directed Logical; should edges be interpreted as directed? Defaults to FALSE.
#' @param loops Logical; should loops be allowed? Defaults to FALSE.
#' @param matrix.type Either "adjacency" or "edgelist".
#' @param vertex.names Optional. A vector or list of vertex names.
#' @param vertex.attr Optional. Additional vertex attributes.
#' @param dual.sign Logical. Allow positive and negative edges simultaneously between the same pair.
#' @param timepoints Optional. Pooling definition for dynamic networks.
#' @param tie.breaker How to resolve ties when pooling signed matrices.
#' @param ... Additional arguments passed to `network::network`.
#'
#' @return A signed network of class `static.sign` or `dynamic.sign`.
#' @export
network.sign <- function(mat = NULL,
                         pos.mat = NULL,
                         neg.mat = NULL,
                         directed = FALSE,
                         loops = FALSE,
                         matrix.type = c("adjacency", "edgelist"),
                         vertex.names = NULL,
                         vertex.attr = NULL,
                         dual.sign = FALSE,
                         timepoints = NULL,
                         tie.breaker = c("zero", "positive", "negative", "first", "last"),
                         ...) {

  matrix.type <- match.arg(matrix.type)
  tie.breaker <- match.arg(tie.breaker)

  convert_to_adj <- function(edgelist, n) {
    adj <- matrix(0, nrow = n, ncol = n)
    if (nrow(edgelist) == 0) return(adj)
    for (i in seq_len(nrow(edgelist))) {
      from <- as.integer(edgelist[i, 1])
      to   <- as.integer(edgelist[i, 2])
      sign <- as.integer(edgelist[i, 3])
      adj[from, to] <- sign
    }
    adj
  }

  pool_timepoints <- function(mats, timepoints, tie.breaker = "zero") {
    stopifnot(all(sapply(mats, function(x) all(x %in% c(-1, 0, 1)))))

    pooled <- vector("list", length(timepoints))

    for (i in seq_along(timepoints)) {
      group <- timepoints[[i]]
      group_mats <- mats[group]

      if (matrix.type == "edgelist") {
        n <- max(sapply(group_mats, function(x) max(as.integer(x[, 1:2]))))
        group_mats <- lapply(group_mats, convert_to_adj, n = n)
      }

      mat_dim <- dim(group_mats[[1]])
      mat_array <- array(unlist(group_mats), dim = c(mat_dim[1], mat_dim[2], length(group_mats)))

      neg_count <- apply(mat_array, 1:2, function(x) sum(x == -1))
      pos_count <- apply(mat_array, 1:2, function(x) sum(x == 1))

      result <- matrix(0, nrow = mat_dim[1], ncol = mat_dim[2])
      result[pos_count > neg_count] <- 1
      result[neg_count > pos_count] <- -1

      tie_mask <- pos_count == neg_count & pos_count > 0
      if (any(tie_mask)) {
        tie_resolved <- matrix(0, nrow = mat_dim[1], ncol = mat_dim[2])
        if (tie.breaker == "zero") {
          tie_resolved[tie_mask] <- 0
        } else if (tie.breaker == "positive") {
          tie_resolved[tie_mask] <- 1
        } else if (tie.breaker == "negative") {
          tie_resolved[tie_mask] <- -1
        } else if (tie.breaker == "first") {
          tie_resolved[tie_mask] <- group_mats[[1]][tie_mask]
        } else if (tie.breaker == "last") {
          tie_resolved[tie_mask] <- group_mats[[length(group_mats)]][tie_mask]
        } else {
          stop("Unknown tie.breaker: ", tie.breaker)
        }
        result[tie_mask] <- tie_resolved[tie_mask]
      }

      ref_dimnames <- dimnames(group_mats[[1]])
      if (!is.null(ref_dimnames[[1]]) && !is.null(ref_dimnames[[2]])) {
        dimnames(result) <- ref_dimnames
      }
      pooled[[i]] <- result
    }
    return(pooled)
  }

  # --- CASE 1: separate positive/negative matrices ---
  if ((!is.null(pos.mat) || !is.null(neg.mat)) && is.null(timepoints)) {
    if (is.null(pos.mat) || is.null(neg.mat))
      stop("Both `pos.mat` and `neg.mat` must be supplied when using dual matrices.")

    if (xor(is.list(pos.mat), is.list(neg.mat)))
      stop("`pos.mat` and `neg.mat` must both be lists or both be single matrices.")

    # --- dynamic ---
    if (is.list(pos.mat)) {
      if (length(pos.mat) != length(neg.mat))
        stop("`pos.mat` and `neg.mat` lists must be of equal length.")

      nets <- mapply(function(p, n) {
        pos_net <- network::network(p, directed = directed, loops = loops,
                                    matrix.type = "adjacency", vertex.attr = vertex.attr, ...)
        neg_net <- network::network(n, directed = directed, loops = loops,
                                    matrix.type = "adjacency", vertex.attr = vertex.attr, ...)
        Layer(`+` = pos_net, `-` = neg_net)
      }, pos.mat, neg.mat, SIMPLIFY = FALSE)

      nets <- lapply(nets, function(net) {
        class(net) <- c("static.sign", class(net))
        net
      })

      MultiDyn <- networks.sign(nets, dynamic = TRUE, dual.sign = dual.sign)
      return(MultiDyn)
    }

    # --- static ---
    pos_net <- network::network(pos.mat, directed = directed, loops = loops,
                                matrix.type = "adjacency", vertex.attr = vertex.attr, ...)
    neg_net <- network::network(neg.mat, directed = directed, loops = loops,
                                matrix.type = "adjacency", vertex.attr = vertex.attr, ...)
    MultiNet <- Layer(`+` = pos_net, `-` = neg_net)

    if (!dual.sign)
      MultiNet %ergmlhs% "constraints" <- update(MultiNet %ergmlhs% "constraints", ~. + fixL(~`+` & `-`))

    MultiNet %v% "sign" <- MultiNet %v% ".LayerName"
    class(MultiNet) <- c("static.sign", "network", class(MultiNet))
    return(MultiNet)
  }

  # --- CASE 1a: dual-sign pooling ---
  if (!is.null(pos.mat) && !is.null(neg.mat) && !is.null(timepoints)) {

    # Ensure lists
    if (!is.list(pos.mat)) pos.mat <- list(pos.mat)
    if (!is.list(neg.mat)) neg.mat <- list(neg.mat)

    if (length(pos.mat) != length(neg.mat))
      stop("`pos.mat` and `neg.mat` lists must be of equal length.")

    if (is.numeric(timepoints)) {
      n_groups <- timepoints
      group_size <- floor(length(pos.mat) / n_groups)
      timepoints <- split(seq_along(pos.mat), ceiling(seq_along(pos.mat) / group_size))
    }

    pooled_pos <- pool_timepoints(pos.mat, timepoints, tie.breaker)
    pooled_neg <- pool_timepoints(neg.mat, timepoints, tie.breaker)

    # Build networks from pooled matrices
    nets <- mapply(function(p, n) {
      pos_net <- network::network(p, directed = directed, loops = loops,
                                  matrix.type = "adjacency", vertex.attr = vertex.attr, ...)
      neg_net <- network::network(n, directed = directed, loops = loops,
                                  matrix.type = "adjacency", vertex.attr = vertex.attr, ...)
      Layer(`+` = pos_net, `-` = neg_net)
    }, pooled_pos, pooled_neg, SIMPLIFY = FALSE)

    nets <- lapply(nets, function(net) {
      class(net) <- c("static.sign", class(net))
      net
    })

    MultiDyn <- networks.sign(nets, dynamic = TRUE, dual.sign = dual.sign)
    return(MultiDyn)
  }

  # --- CASE 2: single signed matrix (original logic) ---
  build_network <- function(x, vertex.names = NULL, vertex.attr = NULL) {
    if (matrix.type == "adjacency") {
      pos_mat <- (x > 0) * 1
      neg_mat <- (x < 0) * 1
      pos_net <- network::network(pos_mat, matrix.type = "adjacency", directed = directed,
                                  loops = loops, vertex.attr = vertex.attr, ...)
      neg_net <- network::network(neg_mat, matrix.type = "adjacency", directed = directed,
                                  loops = loops, vertex.attr = vertex.attr, ...)
      MultiNet <- Layer(`+` = pos_net, `-` = neg_net)

      if (!dual.sign)
        MultiNet %ergmlhs% "constraints" <- update(MultiNet %ergmlhs% "constraints", ~. + fixL(~`+` & `-`))

      MultiNet %v% "sign" <- MultiNet %v% ".LayerName"
      class(MultiNet) <- c("static.sign", "network", class(MultiNet))
      return(MultiNet)

    } else if (matrix.type == "edgelist") {
      x <- as.data.frame(x)
      colnames(x) <- c("from", "to", "sign")
      edges_pos <- x[x[, 3] == 1, ]
      edges_neg <- x[x[, 3] == -1, ]

      pos_net <- network::network(edges_pos, directed = directed, loops = loops,
                                  matrix.type = "edgelist", vertex.attr = vertex.attr, ...)
      neg_net <- network::network(edges_neg, directed = directed, loops = loops,
                                  matrix.type = "edgelist", vertex.attr = vertex.attr, ...)
      MultiNet <- Layer(`+` = pos_net, `-` = neg_net)

      if (!dual.sign)
        MultiNet %ergmlhs% "constraints" <- update(MultiNet %ergmlhs% "constraints", ~. + fixL(~`+` & `-`))

      MultiNet %v% "sign" <- MultiNet %v% ".LayerName"
      class(MultiNet) <- c("static.sign", "network", class(MultiNet))
      return(MultiNet)
    } else {
      stop("Unsupported matrix.type")
    }
  }

  # --- dynamic single-matrix case ---
  if (is.list(mat) && !is.data.frame(mat)) {
    if (!is.null(timepoints)) {
      if (is.numeric(timepoints)) {
        n_groups <- timepoints
        group_size <- floor(length(mat) / n_groups)
        timepoints <- split(seq_along(mat), ceiling(seq_along(mat) / group_size))
      }
      mat <- pool_timepoints(mat, timepoints, tie.breaker)
    }

    if (is.null(vertex.names)) {
      vertex.names <- lapply(mat, function(m) {
        dn <- dimnames(m)
        if (!is.null(dn[[1]])) dn[[1]] else NULL
      })
    } else if (!is.list(vertex.names)) {
      vertex.names <- rep(list(vertex.names), length(mat))
    }

    nets <- mapply(function(x, nms) build_network(x, vertex.names = nms, vertex.attr = vertex.attr),
                   mat,
                   vertex.names,
                   SIMPLIFY = FALSE)

    MultiDyn <- networks.sign(nets, dynamic = TRUE, dual.sign = dual.sign)
    return(MultiDyn)
  } else {
    if (is.null(vertex.names) && matrix.type == "adjacency") {
      dn <- dimnames(mat)
      if (!is.null(dn[[1]])) vertex.names <- dn[[1]]
    }
    return(build_network(mat, vertex.names = vertex.names, vertex.attr = vertex.attr))
  }
}

