#' Fit an ERGM with MPLE using a logistic regression model
#'
#' Return a fitted logistic regression model to calculate the MPLE of an ERGM.
#'
#' The MPLE for an ERGM is calculated by first finding the matrices of
#' positive and negative change statistics. The resulting change matrices are
#' then used to calculate the MPLE using a logistic regression model.
#'
#' @param formula An ERGM formula, with the left-hand side being the network
#' @param ... Additional arguments passed to \code{\link{ergmMPLE}}
#'
#' @return A fitted logistic regression model
#'
#' @seealso \code{\link{ergmMPLE}}
#'
#' @examples
#' data(tribes)
#' mple_sign(tribes ~ Pos(~edges) + Neg(~edges))
#'
#' @export

mple_sign <- function(formula, ...) {
  # Fit the ergmMPLE model
  tmp <- ergmMPLE(formula, output = "array", ...)

  # Extract relevant dimensions
  n_actors <- nrow(tmp$predictor[,,1]) / 2
  n_vars <- dim(tmp$predictor)[3]

  # Extract change statistics matrices
  change_pos <- tmp$predictor[1:n_actors, 1:n_actors, ]
  change_neg <- tmp$predictor[1:n_actors + n_actors, 1:n_actors + n_actors, ]

  # Create adjacency matrix
  net <- ergm.getnetwork(formula)
  adj_mat <- matrix(0, nrow = n_actors, ncol = n_actors)
  adj_mat[which(as.matrix(net)[1:n_actors, 1:n_actors] == 1, arr.ind = TRUE)] <- 1
  adj_mat[which(as.matrix(net)[1:n_actors + n_actors, 1:n_actors + n_actors] == 1, arr.ind = TRUE)] <- -1
  diag(adj_mat) <- NA
  adj_mat[lower.tri(adj_mat)] <- NA

  # Create adjacency array
  adj_array <- array(rep(adj_mat, n_vars), dim = c(n_actors, n_actors, n_vars))

  # Extract indices
  indices_zero <- which(adj_array == 0, arr.ind = TRUE)
  indices_pos <- which(adj_array == 1, arr.ind = TRUE)
  indices_neg <- which(adj_array == -1, arr.ind = TRUE)

  # Extract change statistics based on adjacency values
  zero_1 <- change_pos[indices_zero]
  zero_2 <- change_neg[indices_zero]
  pos_1 <- -change_pos[indices_pos]
  pos_2 <- -change_pos[indices_pos] + change_neg[indices_pos]
  neg_1 <- -change_neg[indices_neg] + change_pos[indices_neg]
  neg_2 <- -change_neg[indices_neg]

  # Convert change statistics to matrices
  mat_zero_1 <- matrix(unlist(zero_1), nrow = ceiling(length(zero_1) / n_vars))
  mat_zero_2 <- matrix(unlist(zero_2), nrow = ceiling(length(zero_2) / n_vars))
  mat_pos_1 <- matrix(unlist(pos_1), nrow = ceiling(length(pos_1) / n_vars))
  mat_pos_2 <- matrix(unlist(pos_2), nrow = ceiling(length(pos_2) / n_vars))
  mat_neg_1 <- matrix(unlist(neg_1), nrow = ceiling(length(neg_1) / n_vars))
  mat_neg_2 <- matrix(unlist(neg_2), nrow = ceiling(length(neg_2) / n_vars))

  # Combine matrices
  result <- rbind(mat_zero_1, mat_zero_2, mat_pos_1, mat_pos_2, mat_neg_1, mat_neg_2)

  # Set column names
  colnames(result) <- dimnames(tmp$predictor)[[3]]

  # Create y variable and bind it to the result matrix
  result <- cbind(y = rep(0, nrow(result)), result)

  # Add additional rows with y = 1 and zeros
  result <- rbind(result, cbind(y = 1, matrix(0, nrow = nrow(result)/2, ncol = ncol(result)-1)))

  # Fit logistic regression model
  fit <- glm(y ~ . - 1, data = as.data.frame(result), family = binomial())

  # Return the fitted model
  return(fit)
}
