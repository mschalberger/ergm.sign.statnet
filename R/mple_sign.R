#' Fit an ERGM with MPLE using a logistic regression model
#'
#' Returns a fitted logistic regression model used to calculate the maximum
#' pseudolikelihood estimate (MPLE) of an exponential random graph model (ERGM).
#'
#' The MPLE is calculated by first computing matrices of positive and negative
#' change statistics. These are then used to estimate the MPLE via logistic
#' regression. Optionally, the covariance can be estimated using the Godambe
#' method.
#'
#' @param formula An ERGM formula with the network on the left-hand side.
#' @param control A list of control parameters for \code{\link[ergm]{ergmMPLE}}.
#'   By default, the covariance method is set to "Godambe".
#' @param seed Optional integer to set the random seed for reproducibility when
#'   simulating networks for Godambe covariance estimation.
#' @param ... Additional arguments passed to \code{\link[ergm]{ergmMPLE}}.
#'
#' @return An object of class \code{\link[ergm]{ergm}}.
#'
#' @seealso \code{\link[ergm]{ergmMPLE}}, \code{\link[ergm]{ergm}}, \code{\link[stats]{glm}}
#'
#' @examples
#' data(tribes)
#' mple_sign(tribes ~ Pos(~edges) + Neg(~edges))
#'
#' @export
mple_sign <- function(formula, control = control.ergm(), seed = NULL, ...) {
  net <- ergm.getnetwork(formula)

  has_fixL <- any(grepl("fixL", deparse(net[["gal"]][["ergm"]][["constraints"]])))

  if (has_fixL) {
    cons_term <- if ("dynamic.sign" %in% class(net)) {
      quote(Cross(~L(~edges, ~`+` & `-`)))
    } else if ("multi.sign" %in% class(net)) {
      quote(N(~L(~edges, ~`+` & `-`)))
    } else {
      quote(L(~edges, ~`+` & `-`))
    }

    formula_use <- update(formula, substitute(. ~ . + TERM, list(TERM = cons_term)))
  } else {
    formula_use <- formula
  }

  # Fit the ergmMPLE model
  mple_data <- ergmMPLE(formula_use, control = control, ...)

  if (has_fixL) {
    keep <- mple_data$predictor[, ncol(mple_data$predictor)] == 0
    predictor <- mple_data$predictor[keep, -ncol(mple_data$predictor), drop = FALSE]
    response  <- mple_data$response[keep]
    weights   <- mple_data$weights[keep]
  } else {
    predictor <- mple_data$predictor
    response  <- mple_data$response
    weights   <- mple_data$weights
  }

  glm_fit <- glm(
    response ~ . - 1,
    data = data.frame(predictor, check.names = FALSE),
    weights = weights,
    family = "binomial"
  )

  glm_summary <- summary(glm_fit)
  res <- list()

  # Godambe covariance estimation
  if (control$MPLE.covariance.method == "Godambe") {
    sim_mple <- ergm::simulate_formula(
      object = formula,
      basis = net,
      nsim = control$MPLE.covariance.samplesize,
      coef = glm_fit$coefficients,
      seed = seed,
      control = control.simulate.formula(
        MCMC.prop = ~sparse,
        MCMC.burnin = control$MPLE.covariance.sim.burnin,
        MCMC.interval = control$MPLE.covariance.sim.interval
      ),
      output = "network",
      ...
    )

    num_variables <- length(glm_fit$coefficients)
    nsim <- length(sim_mple)
    gradient_matrix <- matrix(0, nrow = num_variables, ncol = nsim)

    for (i in seq_len(nsim)) {
      sim_net <- sim_mple[[i]]
      dat <- ergm::ergmMPLE(formula = formula, basis = sim_net, control = control,
                            output = "dyadlist")

      covariates <- dat$predictor[, -(1:2), drop = FALSE]
      predictions <- 1 / (1 + exp(-covariates %*% glm_fit$coefficients))
      gradient <- t(covariates) %*% (dat$response - predictions)
      gradient_matrix[, i] <- gradient
    }

    variability_matrix <- var(t(gradient_matrix))
    invHess <- glm_summary$cov.unscaled
    res$covar <- invHess %*% variability_matrix %*% invHess
  } else {
    res$covar <- glm_summary$cov.unscaled
  }

  glm_fit_null <- glm(
    response ~ 1,
    data = data.frame(predictor, check.names = FALSE),
    weights = weights,
    family = "binomial"
  )



  res$network <- net
  res$coefficients <- glm_fit$coefficients
  res$iterations <- glm_fit$iter
  res$MCMCtheta <- glm_fit$coefficients
  res$gradient <- rep(NA, length(glm_fit$coefficients))
  res$hessian <- -solve(glm_summary$cov.unscaled)
  res$failure <- !glm_fit$converged
  res$mple.lik <- logLik(glm_fit)
  res$mple.lik.null <- logLik(glm_fit_null)
  res$mle.lik <- logLik(glm_fit)
  res$null.lik <- logLik(glm_fit_null)
  res$estimate <- "MPLE"
  res$control <- control
  res$call <- match.call()
  res$constraints <- net[["gal"]][["ergm"]][["constraints"]]
  res$ergm_version <- as.package_version(as.character(packageVersion("ergm")))
  res$formula <- formula
  res$info <- list(
    terms_dind = FALSE,
    space_dind = TRUE,
    n_info_dyads = sum(mple_data$weight),
    obs = FALSE,
    valued = FALSE
  )
  res$etamap$offsettheta <- rep(FALSE, length(res$coefficients))
  res$glm <- glm_fit

  class(res) <- c("ergm")

  return(res)
}

