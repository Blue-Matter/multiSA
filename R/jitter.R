
#' Jitter starting values from fitted model
#'
#' Run additional model fits with jittered starting values.
#'
#' @param x [MSAassess-class] object returned by [fit_MSA()]
#' @param n Integer, number of jittered model runs
#' @param use_fitted Logical, whether to jitter from estimated parameters (`TRUE`) or original starting values (`FALSE`)
#' @param return_models Logical, whether to return fitted models of the jitter runs
#' @param amount Numeric or NULL, passed to [base::jitter()]
#' @param cores Integer, number of CPUs for parallel processing
#' @param seed Integer, for replicating the sampling function. Optional.
#' @param ... Other arguments to pass to [fit_MSA()]
#' @details The new starting parameters are: `pars + amount * runif(n, -1, 1)` where `pars` are either the fitted values or original starting
#' values depending on `use_fitted`.
#' @returns If `return_models = TRUE`, a list (length `n`) containing [MSAassess-class] objects.
#' Otherwise, a data frame of likelihood components made by [get_likelihood_components()]
#' @export
#' @importFrom parallel parLapplyLB
do_jitter <- function(x, n = 1, use_fitted = TRUE, return_models = TRUE, amount = NULL, cores = 1, seed, ...) {

  if (use_fitted) {
    pars <- x@obj$env$last.par.best
  } else {
    pars <- x@obj$par
  }

  if (!missing(seed)) set.seed(seed)
  jit <- lapply(seq_len(n), function(...) jitter(pars, amount = amount))

  jitter_fn <- function(i, fitted, return_models, ...) {

    dat <- get_MSAdata(fitted)
    parameters <- fitted@obj$env$parList(par = i)

    map <- dat@Misc$map
    random <- dat@Misc$random

    fit <- fit_MSA(dat, parameters, map, random, silent = TRUE, ...)

    if (return_models) {
      return(fit)
    } else {
      return(get_likelihood_components(fit))
    }
  }

  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    fits <- parLapplyLB(cl, X = jit, jitter_fn, fitted = x, return_models = return_models, ...)
  } else {
    fits <- lapply(cl, X = jit, jitter_fn, fitted = x, return_models = return_models, ...)
  }

  if (return_models) {
    return(fits)
  } else {
    out <- cbind(
      data.frame(Run = 1:n),
      do.call(rbind, fits)
    )
    return(out)
  }


#
  #init <- fit_MSA(
  #  dat,
  #  parameters = parameters,
  #  map = map,
  #  random = random,
  #  run_model = FALSE
  #)
#
  #init@obj$par[] <- jit
  #fit <- fit_MSA(init, ...)
  #return(fit)
}
