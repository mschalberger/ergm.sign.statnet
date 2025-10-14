#' Create Signed Network Object
#'
#' Turn an adjacency matrix or an edgelist into a static or dynamic signed network.
#'
#' @param mat (List of) adjacency matrices or edgelists. For dynamic networks, provide a list.
#' Adjacency matrices must contain only -1, 0, or 1. Edgelists must have three columns: "From", "To", and "Sign".
#' @param directed Logical; should edges be interpreted as directed? Defaults to FALSE.
#' @param loops Logical; should loops be allowed? Defaults to FALSE.
#' @param matrix.type Either "adjacency" or "edgelist", indicating the format of the input.
#' @param vertex.names Optional. A vector or list of vertex names. If dynamic, provide a list aligned with `mat`.
#' @param vertex.attr Optional. Additional vertex attributes passed to the network constructor.
#' @param timepoints Optional. Either a single integer (number of desired pooled timepoints) or a list of index groups (e.g., list(1:3, 4:6)). Only applies to dynamic networks.
#' @param tie.breaker Character. Determines how to resolve sign conflicts when positives and negatives are equal within a pooling window. One of:
#' \itemize{
#'   \item `"zero"`: set to 0 (default)
#'   \item `"positive"`: set to 1
#'   \item `"negative"`: set to -1
#'   \item `"first"`: use the first nonzero sign
#'   \item `"last"`: use the last nonzero sign
#' }
#' @param ... Additional arguments passed to \code{network::network}.
#'
#' @return A signed network of class \code{"static.sign"} or \code{"dynamic.sign"}.
#'
#' @examples
#' # Static adjacency matrix
#' mat <- matrix(c(0, 1, -1, 1, 0, 0, -1, 0, 0), nrow = 3)
#' net <- signNetwork(mat, matrix.type = "adjacency")
#'
#' # Dynamic network with pooling
#' mats <- replicate(6, mat, simplify = FALSE)
#' net_dyn <- signNetwork(mats, matrix.type = "adjacency", timepoints = 2, tie.breaker = "positive")
#'
#' # Custom timepoint pooling
#' net_custom <- signNetwork(mats, matrix.type = "adjacency", timepoints = list(1:3, 4:6), tie.breaker = "first")
#'
#' @export
signNetwork <- function(mat, directed = FALSE, loops = FALSE, matrix.type = c("adjacency", "edgelist"),
                        vertex.names = NULL, vertex.attr = NULL, timepoints = NULL,
                        tie.breaker = c("zero", "positive", "negative", "first", "last"), ...) {

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

      # Stack into array: [nodes x nodes x timepoints]
      mat_dim <- dim(group_mats[[1]])
      mat_array <- array(unlist(group_mats), dim = c(mat_dim[1], mat_dim[2], length(group_mats)))

      # Count -1s, 0s, and 1s
      neg_count <- apply(mat_array, 1:2, function(x) sum(x == -1))
      pos_count <- apply(mat_array, 1:2, function(x) sum(x == 1))

      result <- matrix(0, nrow = mat_dim[1], ncol = mat_dim[2])

      # Majority decisions
      result[pos_count > neg_count] <- 1
      result[neg_count > pos_count] <- -1

      # Ties (only where pos == neg and > 0)
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

  build_network <- function(x, vertex.names = NULL, vertex.attr = NULL) {
    if (matrix.type == "adjacency") {
      static_net <- network::network(abs(x), matrix.type = "adjacency", directed = directed, loops = loops, vertex.attr = vertex.attr, ...)
      static_net <- network::set.edge.value(static_net, "sign", x)
      static_net%e%"pos" <- static_net%e%"sign" == 1
      static_net%e%"neg" <- static_net%e%"sign" == -1
    } else if (matrix.type == "edgelist") {
      x <- as.data.frame(x)
      colnames(x) <- c("from", "to", "sign")
      edges_pos <- x[x[,3] == 1,]
      edges_neg <- x[x[,3] == -1,]

      edges_pos$pos <- TRUE
      edges_pos$neg <- FALSE
      edges_neg$pos <- FALSE
      edges_neg$neg <- TRUE

      all_edges <- rbind(edges_pos, edges_neg)
      static_net <- network::network(all_edges, matrix.type = "edgelist", directed = directed, loops = loops, vertex.attr = vertex.attr, ...)
    } else {
      stop("Unsupported matrix.type")
    }

    if (!is.null(vertex.names)) {
      network::set.vertex.attribute(static_net, "vertex.names", vertex.names)
    }

    MultiNet <- Layer(static_net, c(`+` = "pos", `-` = "neg"))
    MultiNet %ergmlhs% "constraints" <- update(MultiNet %ergmlhs% "constraints", ~. + fixL(~`+`&`-`))
    MultiNet %v% "sign" <- MultiNet %v% ".LayerName"
    class(MultiNet) <- c("static.sign", "network", class(MultiNet))
    return(MultiNet)
  }

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

    MultiDyn <- signNetworks(nets, dynamic = T)
    #MultiDyn %ergmlhs% "constraints" <- update(MultiDyn %ergmlhs% "constraints", ~ . + fixL(~`+` & `-`))
    MultiDyn$gal$NetList <- nets
    class(MultiDyn) <- c("dynamic.sign", class(MultiDyn))
    return(MultiDyn)
  } else {
    if (is.null(vertex.names) && matrix.type == "adjacency") {
      dn <- dimnames(mat)
      if (!is.null(dn[[1]])) vertex.names <- dn[[1]]
    }
    return(build_network(mat, vertex.names = vertex.names, vertex.attr = vertex.attr))
  }
}

