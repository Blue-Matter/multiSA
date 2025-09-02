

#' @name simulate
#' @aliases simulate.MSAassess
#'
#' @title Simulate data
#'
#' @description Simulate data observations from fitted MSA model.
#'
#' @param object [MSAassess-class] object returned by [fit_MSA()]
#' @param nsim Integer, number of simulations
#' @param seed Random number generator seed
#' @param ... Not used
#' @return A list of `nsim` length with data observations
#' @importFrom stats simulate
#' @export
simulate.MSAassess <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  sims <- lapply(1:nsim, function(...) object@obj$simulate(object@obj$env$last.par.best))

  var_data <- names(object@obj$env$obs)
  out <- lapply(sims, function(x) x[var_data])
  return(out)
}
