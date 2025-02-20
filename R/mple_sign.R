#' Fit an ERGM with MPLE using a logistic regression model
#'
#' Return a fitted logistic regression model to calculate the MPLE of an ERGM.
#'
#' The MPLE for an ERGM is calculated by first finding the matrices of
#' positive and negative change statistics. The resulting change matrices are
#' then used to calculate the MPLE using a logistic regression model.
#'
#' @param formula An ERGM formula, with the left-hand side being the network
#' @param control A list of control parameters for the \code{\link{ergmMPLE}} function. Per default the method to estimate the convariance method is set to Godambe.
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

mple_sign <- function(formula, control = control.ergm(MPLE.covariance.method="Godambe"), ...) {
  # Fit the ergmMPLE model
  tmp <- ergmMPLE(formula, output = "array", control = control, ...)

  net <- ergm.getnetwork(formula)
  n_actors <- network.size(net)/ 2
  n_vars <- dim(tmp$predictor)[3]

  is_directed <- is.directed(net)
  has_loops <- has.loops(net)

  if ("combined_networks" %in% class(net)) {
    k <- length(net[["gal"]][[".subnetcache"]][[".NetworkID"]])
    block_sizes <- unlist(lapply(net[["gal"]][[".subnetcache"]][[".NetworkID"]], function(x) network.size(x)))
    cum_block_sizes <- c(0, cumsum(block_sizes))
    half_sizes <- block_sizes / 2

    change_pos <- array(0, dim = c(sum(half_sizes), sum(half_sizes), n_vars))
    change_neg <- array(0, dim = c(sum(half_sizes), sum(half_sizes), n_vars))
    adj_mat <- matrix(0, nrow = n_actors, ncol = n_actors)


    for (i in 1:k) {
      # Compute the start/end of each positive and negative block
      pos_idx <- (cum_block_sizes[i] + 1):(cum_block_sizes[i] + half_sizes[i])  # First half
      neg_idx <- (cum_block_sizes[i] + half_sizes[i] + 1):cum_block_sizes[i + 1]  # Second half

      # Compute the position in the final extracted matrices
      new_pos_idx <- ((i - 1) * half_sizes[i] + 1):(i * half_sizes[i])

      change_pos[new_pos_idx, new_pos_idx, ] <- tmp$predictor[pos_idx, pos_idx, , drop = FALSE]
      change_neg[new_pos_idx, new_pos_idx, ] <- tmp$predictor[neg_idx, neg_idx, , drop = FALSE]

      # Create adjacency matrix
      adj_mat[new_pos_idx,new_pos_idx][which(as.matrix(net)[pos_idx, pos_idx] == 1, arr.ind = TRUE)] <- 1
      adj_mat[new_pos_idx,new_pos_idx][which(as.matrix(net)[neg_idx, neg_idx] == 1, arr.ind = TRUE)] <- -1
      if(!has_loops) {
        diag(adj_mat) <- NA
        }
      if(!is_directed) {
        adj_mat[lower.tri(adj_mat)] <- NA
        #change_pos[rep(lower.tri(matrix(NA, nrow = dim(change_pos)[1], ncol = dim(change_pos)[2])), dim(change_pos)[3])] <- NA
        #change_neg[rep(lower.tri(matrix(NA, nrow = dim(change_neg)[1], ncol = dim(change_neg)[2])), dim(change_neg)[3])] <- NA
      }
      #adj_mat[new_pos_idx,- new_pos_idx] <- NA
      #adj_mat[-new_pos_idx,new_pos_idx] <- NA
     }
  } else {
  # Extract change statistics matrices
  change_pos <- tmp$predictor[1:n_actors, 1:n_actors, , drop = FALSE]
  change_neg <- tmp$predictor[1:n_actors + n_actors, 1:n_actors + n_actors, , drop = FALSE]

  # Create adjacency matrix
  adj_mat <- matrix(0, nrow = n_actors, ncol = n_actors)
  adj_mat[which(as.matrix(net)[1:n_actors, 1:n_actors] == 1, arr.ind = TRUE)] <- 1
  adj_mat[which(as.matrix(net)[1:n_actors + n_actors, 1:n_actors + n_actors] == 1, arr.ind = TRUE)] <- -1
  if(!has_loops) {
    diag(adj_mat) <- NA
  }
  if(!is_directed) {
    adj_mat[lower.tri(adj_mat)] <- NA
    #change_pos[rep(lower.tri(matrix(NA, nrow = dim(change_pos)[1], ncol = dim(change_pos)[2])), dim(change_pos)[3])] <- NA
    #change_neg[rep(lower.tri(matrix(NA, nrow = dim(change_neg)[1], ncol = dim(change_neg)[2])), dim(change_neg)[3])] <- NA
  }
  }
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

  # Combine matrices by rows (some matrices might be empty)
  matrices <- list(mat_zero_1, mat_zero_2, mat_pos_1, mat_pos_2, mat_neg_1, mat_neg_2)
  list_without_empty_matrices <- Filter(function(x) !(is.matrix(x) && nrow(x) == 0 && ncol(x) == 0), matrices)
  result <- do.call("rbind", list_without_empty_matrices)

  # Set column names
  colnames(result) <- dimnames(tmp$predictor)[[3]]

  # Create y variable and bind it to the result matrix
  result <- cbind(y = rep(0, nrow(result)), result)

  # Add additional rows with y = 1 and zeros
  result <- rbind(result, cbind(y = 1, matrix(0, nrow = nrow(result)/2, ncol = ncol(result)-1)))

  # Fit logistic regression model
  glm_fit <- glm(y ~ . - 1, data = as.data.frame(result), family = binomial())

  glm_fit_null <- glm(y~ 1, family = binomial(), data = as.data.frame(result))
  glm_summary <- summary(glm_fit)
  res <- list()
  res$coefficients <- glm_fit$coefficients
  res$iterations <- glm_fit$iter
  res$MCMCtheta <- glm_fit$coefficients
  res$gradient <- rep(NA,length(glm_fit$coefficients))
  res$hessian <- -solve(glm_summary$cov.unscaled)
  res$covar <- glm_summary$cov.unscaled
  res$failure <- !glm_fit$converged
  res$mple.lik <- logLik(glm_fit)
  res$mple.lik.null <- logLik(glm_fit_null)
  res$mle.lik <- logLik(glm_fit)
  res$null.lik <- logLik(glm_fit_null)
  res$estimate <- "MPLE"
  res$control <- control
  res$ergm_version <- as.package_version("4.6.0")
  res$formula <- formula
  res$info <- list(terms_dind = FALSE, space_dind = TRUE, n_info_dyads = sum(tmp$weight), obs = FALSE,valued = FALSE)
  # res$call <- est$call
  # res$info <- est$info
  res$etamap$offsettheta <- rep(FALSE,length(res$coefficients))
  class(res) <- "ergm"

  # Return the fitted model
  return(res)
}
