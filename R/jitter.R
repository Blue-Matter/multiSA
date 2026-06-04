
#' Jitter starting values from fitted model
#'
#' Run another model fit with jittered starting values. Do multiple jittered fits by calling `jitter()` multiple times, but vary either
#' `amount` or `seed`.
#'
#' @param x [MSAassess-class] object returned by [fit_MSA()]
#' @param use_fitted Logical, whether to jitter from estimated parameters (`TRUE`) or original starting values (`FALSE`)
#' @param amount Numeric or NULL, passed to [base::jitter()]
#' @param seed Integer, for replicating the sampling function. Optional.
#' @param ... Other arguments to pass to [fit_MSA()]
#' @details The new starting parameters are: `pars + amount * runif(n, -1, 1)` where `pars` are either the fitted values or original starting
#' values depending on `use_fitted`.
#' @returns [MSAassess-class] object
#' @export
do_jitter <- function(x, use_fitted = TRUE, amount = NULL, seed, ...) {

  if (use_fitted) {
    pars <- x@obj$env$last.par.best
  } else {
    pars <- x@obj$par
  }

  if (!missing(seed)) set.seed(seed)
  jit <- jitter(pars, amount = amount)

  x@obj$par[] <- jit
  fit <- fit_MSA(x, ...)
  return(fit)
}
