# Function to flatten nested formulas by distributing operations over addition
flatten_formula <- function(f) {
  # Helper function to recursively process expression trees
  flatten_expr <- function(expr) {
    if (length(expr) == 1) {
      # Base case: atomic element (variable name)
      return(list(expr))
    }

    if (is.call(expr) && expr[[1]] == as.name("+")) {
      # If current node is addition, flatten both sides and combine
      left <- flatten_expr(expr[[2]])
      right <- flatten_expr(expr[[3]])
      return(c(left, right))
    } else if (is.call(expr)) {
      # For other operations (like Cross, Pos, gwese, etc.)
      op <- expr[[1]]

      # Find which argument contains a formula (starts with ~)
      formula_arg_idx <- NULL
      for (i in 2:length(expr)) {
        if (is.call(expr[[i]]) && expr[[i]][[1]] == as.name("~")) {
          formula_arg_idx <- i
          break
        }
      }

      if (is.null(formula_arg_idx)) {
        # No formula argument found, return as-is
        return(list(expr))
      }

      # Extract the formula argument and flatten it
      formula_arg <- expr[[formula_arg_idx]][[2]]  # Get RHS of ~
      flattened_args <- flatten_expr(formula_arg)

      # Apply the operation to each flattened argument, preserving other arguments
      return(lapply(flattened_args, function(a) {
        new_call <- expr
        new_call[[formula_arg_idx]] <- call("~", a)
        new_call
      }))
    } else {
      # Not a call, return as-is
      return(list(expr))
    }
  }

  # Check if formula has LHS (e.g., y ~ x)
  has_lhs <- length(f) == 3

  # Extract the right-hand side of the formula
  rhs <- if (has_lhs) f[[3]] else f[[2]]

  # Flatten the expression
  flattened_list <- flatten_expr(rhs)

  # Combine all terms with addition
  result <- flattened_list[[1]]
  if (length(flattened_list) > 1) {
    for (i in 2:length(flattened_list)) {
      result <- call("+", result, flattened_list[[i]])
    }
  }

  # Return as a formula, preserving LHS if present
  if (has_lhs) {
    lhs <- f[[2]]
    as.formula(call("~", lhs, result))
  } else {
    as.formula(call("~", result))
  }
}

path_sampling <- function(net, formula, coef, coef_indep, bridges , nsim = 1, seed = NULL, verbose = 0, control = control.simulate.formula()){
  # Taken from the ergm package
  from = coef_indep
  to = coef
  # Taken from ergm see https://github.com/statnet/ergm/blob/master/R/ergm.bridge.R
  mkpath <- function(n, shift = 0, reverse = FALSE) {
    stopifnot(shift >= -1/2, shift <= 1/2)
    u0 <- seq(from = 0 + 1 / 2 / n, to = 1 - 1 / 2 / n, length.out = n)
    u <- u0 + shift / n
    if (reverse) u <- rev(u)
    list(
      theta = t(rbind(sapply(u, function(u) cbind(to * u + from * (1 - u))))),
      u = u
    )
  }

  path <- mkpath(bridges, shift = 0, reverse = F)
  Dtheta.Du <- (to-from)/bridges
  # Turn the formula into terms and network data
  net <- ergm.getnetwork(formula)
  is_directed <- network::is.directed(net)
  has_loop <- network::has.loops(net)
  # Get the number of actors
  n_actors = network.size(net)/2
  global_stats = summary_formula(formula)
  llrs = numeric(length = bridges)
  llrs = lapply(X = seq_len(bridges), FUN = function(x){
    theta <- path$theta[x, ]

    res <- ergm::simulate_formula(
        object = formula,
        basis = net,
        nsim = nsim,
        coef = theta,
        seed = seed + x,
        control = modifyList(control, list(MCMC.prop = ~sparse,
                                           MCMC.burnin = control$MCMC.burnin,
                                           MCMC.interval = control$MCMC.interval)),
        output = "stats",
        verbose = verbose
      )
    return(res %*%Dtheta.Du)
  })
  return((to-from)%*%global_stats - sum(unlist(llrs)))
}

#' Evaluate Log-Likelihood via Path Sampling
#'
#' This function evaluates the log-likelihood of a fitted signed ergm model using path sampling.
#'
#' @param object A fitted signed ergm model object.
#'
#' @return A logLik object containing the evaluated log-likelihood, degrees of freedom, and number of observations.
#'
#' @export
eval_loglik <- function(object) {
  cat("Evaluating log-likelihood...\n")
  formula <- object$formula
  net <- object$network
  coef <- object$coef

  flattened <- flatten_formula(formula)
  idep <- list_rhs.formula(flattened)[is.dyad.independent(flattened, byterm = T)]
  idep <- Filter(Negate(is.null), idep)
  idep_formula <- as.formula(paste("net ~", paste(idep, collapse = " + ")))

  sub_model <- mple_sign(formula = idep_formula)
  clean_names <- function(x) gsub("`", "", x)

  coef_idep <- coef
  coef_idep[] <- 0
  names_sub <- clean_names(names(coef(sub_model)))
  names_idep <- clean_names(names(coef_idep))

  common <- intersect(names_sub, names_idep)
  coef_idep[common] <- coef(sub_model)[match(common, names_sub)]

  loglk <- path_sampling(formula = formula,
                               coef = coef,
                               coef_indep = coef_idep,
                               net = net,
                               nsim = 1,
                               seed = 123,
                               verbose  = 0,
                               bridges = 20,
                               control = control.simulate.formula(MCMC.burnin = 10000*20,
                                                                  MCMC.interval = 10000,
                               ))
  loglik <- structure(
    loglk + logLik(sub_model),
    df = length(coef),
    nobs = network.size(net)/2,
    class = "logLik")
  return(loglik)
}
