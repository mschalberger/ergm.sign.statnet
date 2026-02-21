#' ERGM Wrapper with Signed Network Support
#'
#' \code{ergm} for signed networks. The \code{mple_sign} function
#' is used when \code{estimate = "MPLE"} or to compute initial values for MCMC.
#' For other networks, behavior is identical to \code{ergm::ergm}.
#'
#' @param formula An \code{\link[stats]{formula}} object specifying the ERGM.
#' @param eval_lik Logical indicating whether to evaluate the likelihood using path sampling.
#' @param control A \code{\link[ergm]{control.ergm}} object.
#' @param ... Additional arguments passed to \code{ergm}.
#'
#' @return An \code{\link[ergm]{ergm}} object.
#'
#' @examples
#' \dontrun{
#' ergm.sign(signed_net ~ edges + triangles) # signed network; uses mple_sign
#' }
#'
#' @export
ergm.sign <- function(formula, estimate = "MLE", eval_lik = FALSE, ..., control = control.ergm()) {
  mple <- mple_sign(formula, eval_lik = FALSE, ..., control = control)
  model <- NULL
  if (estimate != "MPLE") {
  if (is.null(control$init)) control$init <- unname(mple$coefficients)
    model <- ergm::ergm(formula, ..., control = control)
  }
  if (estimate == "MPLE" || (!is.null(model) && isTRUE(model[["MPLE_is_MLE"]]))) model <- mple
  if (eval_lik) model$mle.lik <- model$mple.lik <- eval_loglik(model)
  return(model)
}
