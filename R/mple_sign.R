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

  dual.sign <- ifelse(is.null(net%n%"dual.sign"), TRUE, net%n%"dual.sign")

  if (!dual.sign) {
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

  if (!dual.sign) {
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
    u.data <- matrix(0, nrow = R, ncol = num.variables)
    message("Estimating Godambe Matrix using ", R, " simulated networks.")
    theta.mple <- unname(glm_fit$coefficients)
    invHess    <- glm_summary$cov.unscaled
    i          <- 1
    skipped    <- 0

    for (sim_net in sim_mple) {
      message("  Processing simulated network ", i, " of ", R)
      dat <- ergm::ergmMPLE(formula = formula_use, basis = sim_net)

      if (!dual.sign) {
        tmp_keep <- dat$predictor[, ncol(dat$predictor)] == 0
        X <- dat$predictor[tmp_keep, -ncol(dat$predictor), drop = FALSE]
        y <- dat$response[tmp_keep]
        w <- dat$weights[tmp_keep]
      } else {
        X <- dat$predictor
        y <- dat$response
        w <- dat$weights
      }

      predictions       <- 1 / (1 + exp(-X %*% theta.mple))
      weighted_residuals <- w * (y - predictions)
      gradient          <- t(X) %*% weighted_residuals
      Wvec              <- as.vector(w * predictions * (1 - predictions))
      info              <- t(X) %*% (X * Wvec)

      result <- tryCatch(
        solve(info, gradient),
        error = function(e) {
          message(sprintf("    Skipping simulation %d — singular info matrix: %s", i, conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(result)) {
        u.data[i - skipped, ] <- result
      } else {
        skipped <- skipped + 1
      }

      i <- i + 1
    }

    # Trim unused rows if any draws were skipped
    valid_rows <- R - skipped
    if (skipped > 0) {
      message(sprintf("Godambe: %d of %d simulations skipped due to singularity.", skipped, R))
      u.data <- u.data[seq_len(valid_rows), , drop = FALSE]
    }

    res$covar <- if (valid_rows >= 2) var(u.data) else invHess
  } else {
    res$covar <- glm_summary$cov.unscaled
  }

  glm_fit_null <- glm(
    response ~ 1,
    data = data.frame(predictor, check.names = FALSE),
    weights = weights,
    family = "binomial"
  )

  names(glm_fit$coefficients) <- gsub("`", "", names(coef(glm_fit)))

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
  res$constraints <- net%ergmlhs%"constraints"
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
    res$mle.lik  <- res$mple.lik  <- eval_loglik(res)
  } else {
    ll <- logLik(glm_fit)
    res$mle.lik  <- ll
    res$mple.lik <- ll
  }

  class(res) <- c("ergm")

  return(res)
}

