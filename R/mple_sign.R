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

mple_sign <- function(formula, control = control.ergm(), seed = NULL, ...) {
  new_formula <- update(formula, . ~ . + L(~edges, ~ `+` & `-`))
  # Fit the ergmMPLE model
  tmp <- ergmMPLE(new_formula, #output = "array",
                  control = control, ...)
  keep <- (tmp$predictor[, ncol(tmp$predictor)] == 0)


  # net <- ergm.getnetwork(formula)
  # n_actors <- network.size(net)/ 2
  # n_vars <- dim(tmp$predictor)[3]
  #
  # is_directed <- is.directed(net)
  # has_loops <- has.loops(net)
  #
  # if ("combined_networks" %in% class(net)) {
  #   k <- length(net[["gal"]][[".subnetcache"]][[1]])
  #   block_sizes <- sapply(net[["gal"]][[".subnetcache"]][[1]], network.size)
  #   cum_block_sizes <- c(0, cumsum(block_sizes))
  #   half_sizes <- block_sizes / 2
  #
  #   total_half_size <- sum(half_sizes)
  #
  #   change_pos <- array(NA, dim = c(total_half_size, total_half_size, n_vars))
  #   change_neg <- array(NA, dim = c(total_half_size, total_half_size, n_vars))
  #   adj_mat <- matrix(0, nrow = n_actors, ncol = n_actors)
  #
  #   net_mat <- as.matrix(net)  # Convert network to matrix once
  #
  #   new_pos_start <- 1  # Track where the next block starts in the extracted matrices
  #
  #   for (i in seq_len(k)) {
  #     half_size <- half_sizes[i]
  #     pos_idx <- seq(cum_block_sizes[i] + 1, cum_block_sizes[i] + half_size)
  #     neg_idx <- seq(cum_block_sizes[i] + half_size + 1, cum_block_sizes[i + 1])
  #
  #     new_pos_idx <- seq(new_pos_start, new_pos_start + half_size - 1)
  #     new_pos_start <- new_pos_start + half_size  # Update for next iteration
  #
  #     change_pos[new_pos_idx, new_pos_idx, ] <- tmp$predictor[pos_idx, pos_idx, , drop = FALSE]
  #     change_neg[new_pos_idx, new_pos_idx, ] <- tmp$predictor[neg_idx, neg_idx, , drop = FALSE]
  #
  #     adj_mat[new_pos_idx, new_pos_idx][net_mat[pos_idx, pos_idx] == 1] <- 1
  #     adj_mat[new_pos_idx, new_pos_idx][net_mat[neg_idx, neg_idx] == 1] <- -1
  #
  #     adj_mat[new_pos_idx,- new_pos_idx] <- NA
  #     adj_mat[- new_pos_idx, new_pos_idx] <- NA
  #   }
  #
  #   if (!has_loops) {
  #     diag(adj_mat) <- NA
  #   }
  #   if (!is_directed) {
  #     adj_mat[lower.tri(adj_mat)] <- NA
  #   }
  # }
  # else {
  # # Extract change statistics matrices
  # change_pos <- tmp$predictor[1:n_actors, 1:n_actors, , drop = FALSE]
  # change_neg <- tmp$predictor[1:n_actors + n_actors, 1:n_actors + n_actors, , drop = FALSE]
  #
  # # Create adjacency matrix
  # adj_mat <- matrix(0, nrow = n_actors, ncol = n_actors)
  # adj_mat[which(as.matrix(net)[1:n_actors, 1:n_actors] == 1, arr.ind = TRUE)] <- 1
  # adj_mat[which(as.matrix(net)[1:n_actors + n_actors, 1:n_actors + n_actors] == 1, arr.ind = TRUE)] <- -1
  # if(!has_loops) {
  #   diag(adj_mat) <- NA
  # }
  # if(!is_directed) {
  #   adj_mat[lower.tri(adj_mat)] <- NA
  #   #change_pos[rep(lower.tri(matrix(NA, nrow = dim(change_pos)[1], ncol = dim(change_pos)[2])), dim(change_pos)[3])] <- NA
  #   #change_neg[rep(lower.tri(matrix(NA, nrow = dim(change_neg)[1], ncol = dim(change_neg)[2])), dim(change_neg)[3])] <- NA
  # }
  # }
  # # Create adjacency array
  # adj_array <- array(rep(adj_mat, n_vars), dim = c(n_actors, n_actors, n_vars))
  #
  # # Extract indices
  # indices_zero <- which(adj_array == 0, arr.ind = TRUE)
  # indices_pos <- which(adj_array == 1, arr.ind = TRUE)
  # indices_neg <- which(adj_array == -1, arr.ind = TRUE)
  #
  # # Extract change statistics based on adjacency values
  # zero_1 <- change_pos[indices_zero]
  # zero_2 <- change_neg[indices_zero]
  # pos_1 <- -change_pos[indices_pos]
  # pos_2 <- -change_pos[indices_pos] + change_neg[indices_pos]
  # neg_1 <- -change_neg[indices_neg] + change_pos[indices_neg]
  # neg_2 <- -change_neg[indices_neg]
  #
  # # Convert change statistics to matrices
  # mat_zero_1 <- matrix(unlist(zero_1), nrow = ceiling(length(zero_1) / n_vars))
  # mat_zero_2 <- matrix(unlist(zero_2), nrow = ceiling(length(zero_2) / n_vars))
  # mat_pos_1 <- matrix(unlist(pos_1), nrow = ceiling(length(pos_1) / n_vars))
  # mat_pos_2 <- matrix(unlist(pos_2), nrow = ceiling(length(pos_2) / n_vars))
  # mat_neg_1 <- matrix(unlist(neg_1), nrow = ceiling(length(neg_1) / n_vars))
  # mat_neg_2 <- matrix(unlist(neg_2), nrow = ceiling(length(neg_2) / n_vars))
  #
  # # Combine matrices by rows (some matrices might be empty)
  # matrices <- list(mat_zero_1, mat_zero_2, mat_pos_1, mat_pos_2, mat_neg_1, mat_neg_2)
  # list_without_empty_matrices <- Filter(function(x) !(is.matrix(x) && nrow(x) == 0 && ncol(x) == 0), matrices)
  # result <- do.call("rbind", list_without_empty_matrices)
  #
  # # Set column names
  # colnames(result) <- dimnames(tmp$predictor)[[3]]
  #
  # # Create y variable and bind it to the result matrix
  # result <- cbind(y = rep(0, nrow(result)), result)
  #
  # # Add additional rows with y = 1 and zeros
  # result <- rbind(result, cbind(y = 1, matrix(0, nrow = nrow(result)/2, ncol = ncol(result)-1)))

  # Fit logistic regression model
  #glm_fit <- glm(y ~ . - 1, data = as.data.frame(result), family = binomial())
  glm_fit <- glm(tmp$response[keep]  ~ . - 1,
                 data = data.frame(tmp$predictor[keep,-ncol(tmp$predictor),drop = FALSE], check.names = FALSE),
                 weights = tmp$weights[keep],
                 family="binomial")
  glm_summary <- summary(glm_fit)
  res <- list()



  #Godambe
  if (control$MPLE.covariance.method == 'Godambe') {
  if ("combined_networks" %in% class(net)) {
    list_net <- uncombine_network(net)
    small <- signNetwork(matrix(0, nrow = 1, ncol = 1), matrix.type= "adjacency", directed = is_directed, loops = has_loops)
    vnames <- net%v%"vertex.names"
    sub_list <- list()


    if (!is.null(seed)) {
      set.seed(seed)
      seeds <- sample.int(1e6, length(list_net))
    }

    start_name <- 1  # Initialize `start_name`

    for (i in seq_along(list_net)) {
      blocks <- list_net[[i]]
      sample_seed <- if (!is.null(seed)) seeds[[i]] else NULL
      block_size <- network.size(blocks)

      new_names <- vnames[start_name:(start_name + block_size - 1)]
      start_name <- start_name + block_size  # Update correctly

      class(blocks) <- c("static.sign", class(blocks))
      net_list <- signNetworks(blocks, small)

      sub_net <- ergm::simulate_formula(
          object = formula,
          basis = net_list,
          nsim = control$MPLE.covariance.samplesize,  # 500 simulations
          coef = glm_fit$coefficients,
          seed = sample_seed,
          control = control.simulate.formula(
            MCMC.prop = ~sparse,
            MCMC.burnin = control$MPLE.covariance.sim.burnin,
            MCMC.interval = control$MPLE.covariance.sim.interval
          ),
          output = "network",
          ...
      )
      sub_net <- lapply(sub_net, function(net) {
        class(net) <- c("static.sign", class(net))
        net <- get.inducedSubgraph(net, which(net%v%".NetworkID" == 1))
        })
      sub_list[[i]] <- sub_net
    }
    transposed_list <- lapply(seq_len(control$MPLE.covariance.samplesize), function(i) {
      lapply(sub_list, `[[`, i)
    })

    sim_mple <- lapply(transposed_list, signNetworks)
  } else {
    sim_mple <- ergm::simulate_formula(
      object = formula,
      basis = net,
      nsim = control$MPLE.covariance.samplesize,  # 500 simulations
      coef = glm_fit$coefficients,
      control = control.simulate.formula(
        MCMC.prop = ~sparse,
        MCMC.burnin = control$MPLE.covariance.sim.burnin,
        MCMC.interval = control$MPLE.covariance.sim.interval
      ),
      output = "network",
      ...
    )
  }

    num_variables <- length(glm_fit$coefficients)
    nsim <- length(sim_mple)

    gradient_matrix <- matrix(0, nrow = num_variables, ncol = nsim)

    for (i in 1:nsim) {
      tmp <- sim_mple[[i]]
      tmp_match <- tmp %v% "vertex.names"

      dat <- ergm::ergmMPLE(formula = formula, basis = tmp, control = control,
                            output = "dyadlist")
      # Reorder
      dat$predictor[,1] <- tmp_match[dat$predictor[,1]]
      dat$predictor[,2] <- tmp_match[dat$predictor[,2]]

      if(ncol(dat$predictor) == 3){
        covariates <- matrix(ncol = 1, dat$predictor[,-(1:2)])
      } else {
        covariates <- matrix(ncol = ncol(dat$predictor[,-(1:2)]), dat$predictor[,-(1:2)])
      }
      colnames(covariates) <- colnames(dat$predictor)[-(1:2)]

      predictions <- as.vector(1 / (1 + exp(-covariates %*% glm_fit$coef)))

      gradient <- as.vector((dat$response - predictions) %*% covariates)

      gradient_matrix[, i] <- gradient
    }
    # Compute Variability Matrix
    variability_matrix <- var(t(gradient_matrix))

    invHess <- glm_summary$cov.unscaled

    # Godambe matrix calculation
    G <- invHess %*% variability_matrix %*% invHess
    res$covar <- G
  } else {
    res$covar <- glm_summary$cov.unscaled
  }

  glm_fit_null <- glm(tmp$response[keep]  ~  1,
                      data = data.frame(tmp$predictor[keep,-ncol(tmp$predictor), drop = F]),
                      weights = tmp$weights[keep],
                      family="binomial")
  res$coefficients <- glm_fit$coefficients
  res$iterations <- glm_fit$iter
  res$MCMCtheta <- glm_fit$coefficients
  res$gradient <- rep(NA,length(glm_fit$coefficients))
  res$hessian <- -solve(glm_summary$cov.unscaled)
  res$failure <- !glm_fit$converged
  res$mple.lik <- logLik(glm_fit)
  res$mple.lik.null <- logLik(glm_fit_null)
  res$mle.lik <- logLik(glm_fit)
  res$null.lik <- logLik(glm_fit_null)
  res$estimate <- "MPLE"
  res$control <- control
  res$ergm_version <- as.package_version("4.9")
  res$formula <- formula
  res$info <- list(terms_dind = FALSE, space_dind = TRUE, n_info_dyads = sum(tmp$weight), obs = FALSE,valued = FALSE)
  # res$call <- est$call
  # res$info <- est$info
  res$etamap$offsettheta <- rep(FALSE,length(res$coefficients))
  res$glm <- glm_fit
  class(res) <- "ergm"

  # Return the fitted model
  return(res)
}
