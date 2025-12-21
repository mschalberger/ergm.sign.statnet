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
#' @param eval_lik Logical indicating whether to evaluate the likelihood using path sampling.
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
mple_sign <- function(formula, control = control.ergm(), seed = NULL, eval_lik = FALSE,...) {
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
    R <- control$MPLE.covariance.samplesize
    mple.burnin <- control$MPLE.covariance.sim.burnin
    mple.interval <- control$MPLE.covariance.sim.interval

    sim_mple <- ergm::simulate_formula(
      object = formula,
      basis = net,
      nsim = R,
      coef = unname(glm_fit$coefficients),
      control = control.simulate.formula(
        MCMC.prop = ~sparse,
        MCMC.burnin = mple.burnin,
        MCMC.interval = mple.interval
      ),
      output = "network",
      ...
    )

    message("Simulating networks for Godambe covariance estimation...")

    num.variables <- length(glm_fit$coefficients)
    #U <- matrix(0, nrow = R, ncol = num_variables)
    u.data <- matrix(0, nrow = R, ncol = num.variables)
    old.data <- matrix(0, nrow = R, ncol = num.variables)
    message("Estimating Godambe Matrix using ", R, " simulated networks.")
    theta.mple <- unname(glm_fit$coefficients)
    invHess <- glm_summary$cov.unscaled
    i <- 1
    for (sim_net  in sim_mple) {
      message("  Processing simulated network ", i, " of ", R)
      dat <- ergm::ergmMPLE(formula = formula_use, basis = sim_net)
      sim_mple[[i]] <- NULL
      gc(FALSE)

       if (has_fixL) {
         tmp_keep <- dat$predictor[, ncol(dat$predictor)] == 0
         X <- dat$predictor[tmp_keep, -ncol(dat$predictor), drop = FALSE]
         y  <- dat$response[tmp_keep]
         w   <- dat$weights[tmp_keep]
       } else {
         X <- dat$predictor
         y  <- dat$response
         w   <- dat$weights
       }

      rm(dat)
      gc(FALSE)
      predictions <- 1 / (1 + exp(-X %*% theta.mple))
      weighted_residuals <- w * (y - predictions)
      gradient <- t(X) %*% weighted_residuals
      Wvec <- as.vector(w * predictions * (1 - predictions))
      info <- t(X) %*% (X * Wvec)
      u.data[i,] <- solve(info, gradient)
      old.data[i,] <- gradient
      i <- i +1
    }
    res$old_covar <- invHess %*% var(old.data) %*% invHess
    res$covar <- var(u.data)
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
  res$mple.lik.null <- logLik(glm_fit_null)
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
  if (eval_lik) {
    res$mle.lik  <- eval_loglik(glm_fit)
    res$mple.lik <- eval_loglik(glm_fit)
  } else {
    ll <- logLik(glm_fit)
    res$mle.lik  <- ll
    res$mple.lik <- ll
  }

  class(res) <- c("ergm")

  return(res)
}

