#' Create Signed Network Object
#'
#' Turn an adjacency matrix or an edgelist into a static or dynamic signed network
#'
#' @param mat (List of) Adjacency matrix or edgelist. The adjacency matrix must only consist of 1,0 or -1. The edgelist must consist of 3 columns "From", "To" and "Sign" (1 or -1). For a dynamic network the required input is a list of adjacecny matrices or edgelists.
#' @param directed 	logical; should edges be interpreted as directed?
#' @param loops logical; should loops be allowed?
#' @param matrix.type Either "adjacency" or "edgelist" indicating what kind of format the input has.
#' @param cov Add vertex attributes to the network. The input for this should be a dataframe where the first column contains the names of the vertices and the other columns of the dataframe should represent the individual vertex attributes that will be added to the network, however it is not necessary to include all the vertices of the network in the dataframe.
#' @param names Specify vertex names. If no list is provided, the column names of the adjacency matrix will be used as names for the vertices.
#'
#' @return A signed network of the class \code{static.sign} or \code{dynamic.sign}.
#'
#' @export


signnet <- function(mat, directed = F, loops = F, matrix.type, cov = NULL, names = NULL, ...) {
  if (is.list(mat)) {
    i <- 0
    nets <- lapply(mat, function (x) {
      i <- i + 1
      if (matrix.type == "adjacency") {
        if (directed == F){
          if (isSymmetric(x)) {
            x
          } else if (all.equal(x[lower.tri(x)], rep(0,(nrow(x)^2-nrow(x))/2)) == T || #lower triangle is all zeros or
              all.equal(x[lower.tri(x)], rep(NA,(nrow(x)^2-nrow(x))/2)) == T) { #lower triangle is all NAs
            x[is.na(x)] <- 0
            x <- x + t(x)
          } else {
            stop("Matrix is not symmetric!")
          }
        }
        if (loops == F){
          diag(x) <-0
        }
        dyn_net <- as.network(x, matrix.type = "adjacency", directed = directed, loops = loops, ... = ...)
        dyn_net <- set.edge.value(dyn_net, "sign", x)
        pos <- (dyn_net%e%"sign" == 1)
        neg <- (dyn_net%e%"sign" == -1)
        dyn_net%e%"pos" <- pos
        dyn_net%e%"neg" <- neg

      } else if (matrix.type == "edgelist"){
        x <- as.data.frame(x)
        colnames(x) <- c("from","to","sign")
        edges_pos <- x[x[,3]== 1,]
        tryCatch({
          edges_pos$pos = TRUE
          edges_pos$neg = FALSE}, error = function(e){})
        edges_neg <- x[x[,3]== -1,]
        tryCatch({
          edges_neg$pos = FALSE
          edges_neg$neg = TRUE}, error = function(e){})
        full <- rbind(edges_pos, edges_neg)
        dyn_net <- as.network(full, matrix.type = "edgelist", directed = directed, loops = loops, ... = ...)
      }
      if(is.null(names)) {
        if (matrix.type == "adjacency") {
          tryCatch({network.vertex.names(dyn_net) <- rownames(x)}, error = function(e){})
        }
      } else {
        network.vertex.names(dyn_net) <- names[[i]]
      }

      if (!is.null(cov)) {
        # match the names in the network.vertex.names(net) list with the names in the first column of the cov data frame
        match_indices <- match(network.vertex.names(dyn_net), cov[,1])

        # subset the cov data frame based on the match indices
        subset_cov <- cov[match_indices, ]

        # sort the subset_cov data frame based on the network.vertex.names(net) list
        subset_cov <- subset_cov[order(match(subset_cov[,1], network.vertex.names(dyn_net))), ]

        for(col in colnames(subset_cov)[-1]){
          set.vertex.attribute(dyn_net, col, subset_cov[,col])
        }
      }
      x <- dyn_net
    })
    #nets <- networkDynamic(network.list = nets, create.TEAs = T)
    class(nets) <- "dynamic.sign"
    return(nets)
  }
  else if (is.matrix(mat)) {
  if (matrix.type == "adjacency") {
    if (directed == F){
      if (isSymmetric(mat)) {
        mat
      } else if (all.equal(mat[lower.tri(mat)], rep(0,(nrow(mat)^2-nrow(mat))/2)) == T || #lower triangle is all zeros or
          all.equal(mat[lower.tri(mat)], rep(NA,(nrow(mat)^2-nrow(mat))/2)) == T) { #lower triangle is all NAs
        mat[is.na(mat)] <- 0
        mat <- mat + t(mat)
      } else {
        stop("Matrix is not symmetric!")
      }
    }
    if (loops == F){
      diag(mat) <-0
    }

    static_net <- as.network(abs(mat), matrix.type = "adjacency", directed = directed, loops = loops, ... = ...)
    static_net <- set.edge.value(static_net, "sign", mat)
    pos <- (static_net%e%"sign" == 1)
    neg <- (static_net%e%"sign" == -1)
    static_net%e%"pos" <- pos
    static_net%e%"neg" <- neg

  } else if (matrix.type == "edgelist"){
    mat <- as.data.frame(mat)
    colnames(mat) <- c("from","to","sign")
    edges_pos <- mat[mat[,3]== 1,]
    tryCatch({
      edges_pos$pos = TRUE
      edges_pos$neg = FALSE}, error = function(e){})
    edges_neg <- mat[mat[,3]== -1,]
    tryCatch({
      edges_neg$pos = FALSE
      edges_neg$neg = TRUE}, error = function(e){})
    static_net <- as.network(rbind(edges_pos, edges_neg), matrix.type = "edgelist", directed = directed, loops = loops, ... = ...)
  }

  if (is.null(names)) {
    if (matrix.type == "adjacency") {
      tryCatch({network.vertex.names(static_net) <- rownames(mat)}, error = function(e){})
    }
  } else {
    network.vertex.names(static_net) <- names
  }

  if (!is.null(cov)) {
    # match the names in the network.vertex.names(net) list with the names in the first column of the cov data frame
    match_indices <- match(network.vertex.names(static_net), cov[,1])

    # subset the cov data frame based on the match indices
    subset_cov <- cov[match_indices, ]

    # sort the subset_cov data frame based on the network.vertex.names(net) list
    subset_cov <- subset_cov[order(match(subset_cov[,1], network.vertex.names(static_net))), ]

    for(col in colnames(subset_cov)[-1]){
      set.vertex.attribute(static_net, col, subset_cov[,col])
    }
  }

  class(static_net) <- "static.sign"
  return(static_net)
  }
}

