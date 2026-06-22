
#' Quadratic penalty function
#'
#' Taped penalty function if `x < eps`
#'
#' @param x Numeric, the parameter
#' @param eps Numeric, the threshold below which a penalty will be applied
#'
#' @return
#' The penalty value is
#'
#' \deqn{
#' \textrm{penalty} =
#' \begin{cases}
#' 0.1 (x - \varepsilon)^2 & x \le \varepsilon\\
#' 0 & x > \varepsilon
#' \end{cases}
#' }
#'
#' @return Numeric
#' @export
posfun <- function(x, eps) CondExpGe(x, eps, 0, 0.01 * (x - eps) * (x - eps))

calc_selpar_penalty <- function(sel_pf, sel_f, lmid, na, map) {
  if (missing(map) || is.null(map)) {
    map <- array(1:length(sel_pf), dim(sel_pf))
  } else {
    map <- array(map, dim(sel_pf))
  }
  map <- map[-1, ]
  map_num <- as.numeric(map)

  unique_par <- !is.na(map_num) & !duplicated(map_num)

  if (inherits(sel_pf, "advector")) `[<-` <- RTMB::ADoverload("[<-")

  if (all(is.na(map))) {
    if (inherits(sel_pf, "advector")) val <- advector(0) else val <- 0
  } else {
    penalty <- array(0, c(2, ncol(sel_pf)))
    parametric_sel <- grepl("dome|logistic", sel_f)
    flen <- parametric_sel & grepl("length", sel_f)
    if (any(flen)) {
      log_binwidth <- log(0.5 * min(diff(lmid)))
      log_binrange <- log(max(lmid) - min(lmid))
      pen_len <- posfun(sel_pf[2:3, flen], log_binwidth) + posfun(log_binrange, sel_pf[2:3, flen])
      penalty[, flen] <- penalty[, flen] + pen_len
    }
    fage <- parametric_sel & grepl("age", sel_f)
    if (any(fage)) {
      pen_age <- posfun(sel_pf[2:3, fage], log(0.5)) + posfun(log(na), sel_pf[2:3, fage])
      penalty[, fage] <- penalty[, fage] + pen_age
    }
    val <- sum(penalty[unique_par])
  }

  return(val)
}

#' Softmax function
#'
#' Takes a vector of real numbers and returns the corresponding vector of probabilities
#'
#' @param eta Vector
#' @param log Logical, whether to return the value of the logarithm
#'
#' @return Numeric, vector length of `eta`: \eqn{\exp(\eta)/\sum\exp(\eta)}
#' @details Uses `multiSA:::logspace.add` for numerical stability
#' @export
softmax <- function(eta, log = FALSE) {
  den <- Reduce(logspace.add, eta)
  v <- eta - den

  if (log) {
    v
  } else {
    exp(v)
  }
}

logspace.add <- function(lx, ly) CondExpGt(lx, ly, lx, ly) + log1p(exp(-abs(lx - ly)))

#' Calculate covariance matrix
#'
#' Uses Cholesky factorization to generate a covariance matrix (or any symmetric positive definite matrix).
#'
#' @param sigma Numeric vector of marginal standard deviations (all greater than zeros)
#' @param lower_diag Numeric vector to populate the lower triangle of the correlation matrix. All real numbers.
#' Length `sum(1:(length(sigma) - 1))`
#' @examples
#' set.seed(23)
#' n <- 5
#' sigma <- runif(n, 0, 2)
#' lower_diag <- runif(sum(1:(n-1)), -10, 10)
#' Sigma <- conv_Sigma(sigma, lower_diag)
#' Sigma/t(Sigma) # Is symmetric matrix? All ones
#' cov2cor(Sigma)
#'
#' @return Numeric
#' @export
conv_Sigma <- function(sigma, lower_diag) {
  n <- length(sigma)
  stopifnot(length(lower_diag) == sum(1:(n-1)))

  # Parameterizes correlation matrix of X in terms of Cholesky factors
  # https://github.com/kaskr/RTMB/blob/master/tmb_examples/sdv_multi.R
  L <- diag(n)
  L[lower.tri(L)] <- lower_diag
  row_norms <- apply(L, 1, function(row) sqrt(sum(row * row)))
  L <- L / row_norms
  R <- L %*% t(L)  # Correlation matrix of X (guaranteed positive definite)

  V <- diag(sigma)
  Sigma <- V %*% R %*% V
  return(Sigma)
}

conv_steepness <- function(x, SRR = c("BH", "Ricker")) {
  SRR <- match.arg(SRR)
  switch(
    SRR,
    "BH" = 0.8 * plogis(x),
    "Ricker" = exp(x)
  ) + 0.2
}

conv_mat <- function(x, na) {
  a50 <- na * plogis(x[1])
  a95 <- a50 + exp(x[2])

  a <- seq(1, na) - 1
  m <- 1/(1 + exp(-log(19) * (a - a50)/(a95 - a50)))
  return(m)
}


#' Optimize RTMB model
#'
#' A convenient function to fit a RTMB model with [stats::nlminb()]
#'
#' @param obj The list returned by [RTMB::MakeADFun()]
#' @param hessian Logical, whether to pass the Hessian function `obj$he` to [stats::nlminb()]. Only used if
#' there are no random effects in the model.
#' @param restart Deprecated.
#' @param do_sd Deprecated.
#' @param control List of options passed to [stats::nlminb()]
#' @param lower Lower bounds of parameters passed to [stats::nlminb()]
#' @param upper Upper bounds of parameters passed to [stats::nlminb()]
#' @param silent Logical, whether to report progress to console
#' @return A named list, output of [stats::nlminb()]
#' @importFrom stats nlminb
#' @seealso [get_sdreport()]
#' @keywords internal
#' @export
optimize_RTMB <- function(obj, hessian = FALSE, restart = 0, do_sd = TRUE,
                          control = list(iter.max = 2e+05, eval.max = 4e+05),
                          lower = -Inf, upper = Inf, silent = FALSE) {

  if (is.null(obj$env$random) && hessian) {
    h <- obj$he
    if (!silent) message("Using hessian in optimization (can be memory-intensive)")
  } else {
    h <- NULL
  }

  if (!silent) message_info("Fitting model with stats::nlminb()..")
  opt <- tryCatch(
    nlminb(obj$par, obj$fn, obj$gr, h, control = control, lower = lower, upper = upper),
    error = function(e) as.character(e)
  )
  if (!silent) message_info("Final gradient is ", round(max(abs(obj$gr(obj$env$last.par.best))), 5))

  return(opt)
}

# Check that hessian is positive-definite
check_h <- function(h) {
  L <- try(chol(h), silent = TRUE)
  !is.character(L)
}

# Check that hessian could be marginally positive-definite: abs(det(h)) < tol
marginal_h <- function(h, tol = 0.1) {
  det_h <- determinant(h)
  !is.na(det_h$modulus) && det_h$modulus < log(tol)
}

#' Calculate standard errors
#'
#' A wrapper function to return standard errors. Various numerical techniques are employed to obtain
#' a positive-definite covariance matrix in marginal cases.
#' @inheritParams optimize_RTMB
#' @param par.fixed Numeric vector of parameters from which to calculate covariance matrix. Optional
#' @param exact Logical, whether to use autodiff or finite-difference approximation for the hessian. See details.
#' @param getReportCovariance Logical, passed to [RTMB::sdreport()]
#' @param silent Logical, whether to report progress to console. See details.
#' @param ... Other arguments to [RTMB::sdreport()] besides `par.fixed, hessian.fixed, getReportCovariance`
#' @details
#' Uses [stats::optimHess()] if `exact = FALSE`.
#' Autodiff with `exact = TRUE` is only available for TMB models without random effects, but is also memory-intensive.
#'
#' In numerically marginal cases where the determinant of the Hessian matrix is less than 0.1, the function will attempt
#' to calculate the hessian with `numDeriv::jacobian()` and the gradient from TMB.
#'
#' Finally, in other marginal cases where [chol()] identifies a positive-definite Hessian but [solve()] fails to
#' invert the matrix, the covariance matrix will be updated with `chol2inv(chol(h))`
#' @return
#' Object returned by [RTMB::sdreport()].
#'
#' A correlation matrix is generated and stored in: `get_sdreport(obj)$env$corr.fixed`
#'
#' The hessian is stored in `get_sdreport(obj)$env$hessian`
#' @importFrom stats optimHess
#' @export
get_sdreport <- function(obj, par.fixed, exact = FALSE, getReportCovariance = FALSE, silent = FALSE, ...) {
  if (missing(par.fixed)) {
    par.fixed <- obj$env$last.par.best
    if (!is.null(obj$env$random)) par.fixed <- par.fixed[-obj$env$random]
  }

  if (exact && is.null(obj$env$random)) {
    if (!silent) message_info("Calculating standard errors with hessian from obj$he()..")
    h <- obj$he(par.fixed)
  } else {
    if (!silent) message_info("Calculating standard errors with hessian from stats::optimHess()..")
    h <- optimHess(par.fixed, obj$fn, obj$gr)
  }
  res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                  getReportCovariance = getReportCovariance, ...)

  if (!res$pdHess) {
    if (!silent) message_oops("Hessian is not positive-definite.")

    if (marginal_h(h) && requireNamespace("numDeriv", quietly = TRUE)) {
      if (!silent) message_info("Calculating standard errors with hessian from numDeriv::jacobian()..")
      h <- numDeriv::jacobian(obj$gr, par.fixed)
      h <- 0.5 * (h + t(h)) # glmmTMB does this

      if (check_h(h)) {
        res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                        getReportCovariance = getReportCovariance, ...)
      } else if (!silent) {
        message_oops("Hessian is not positive-definite.")
      }
    }
  }

  if (any(is.na(res$cov.fixed)) && res$pdHess) {
    if (!silent) message_info("Calculating standard errors from chol2inv(chol(h))..")
    ch <- try(chol(h), silent = TRUE) # Not needed, this is the test for convergence in sdreport
    if (!is.character(ch)) res$cov.fixed <- chol2inv(ch)
  }

  fixed.names <- make_unique_names(res, select = "fixed")

  res$env$corr.fixed <- cov2cor(res$cov.fixed) |> round(3) |>
    structure(dimnames = list(fixed.names, fixed.names))

  res$env$hessian <- round(h, 3) |>
    structure(dimnames = list(fixed.names, fixed.names))

  return(res)
}

#' @importFrom TMB summary.sdreport
sdreport_int <- function(object, select = c("all", "fixed", "random", "report"), p.value = FALSE, ...) {
  if (is.character(object)) return(object)
  select <- match.arg(select, several.ok = TRUE)
  if ("all" %in% select) select <- c("fixed", "random", "report")
  if ("report" %in% select) {
    AD <- TMB::summary.sdreport(object, "report", p.value = p.value) |> cbind("Gradient" = NA_real_)
    ADnames <- make_unique_names(object, select = "report")
  } else AD <- ADnames <- NULL

  if ("fixed" %in% select) {
    fix <- TMB::summary.sdreport(object, "fixed", p.value = p.value) |> cbind("Gradient" = as.vector(object$gradient.fixed))
    fixnames <- make_unique_names(object, select = "fixed")
  } else fix <- fixnames <- NULL

  if (!is.null(object$par.random) && "random" %in% select) {
    random <- TMB::summary.sdreport(object, "random", p.value = p.value) |> cbind("Gradient" = rep(NA_real_, length(object$par.random)))
    randomnames <- make_unique_names(object, select = "random")
  } else {
    random <- randomnames <- NULL
  }

  out <- rbind(AD, fix, random)
  out <- cbind(out, "CV" = ifelse(abs(out[, "Estimate"]) > 0, out[, "Std. Error"]/abs(out[, "Estimate"]), NA_real_))
  rownames(out) <- c(ADnames, fixnames, randomnames)
  return(out)
}


#' Retrieve data object used to fit model
#'
#' A convenient function to retrieve the data object used to fit the model. The object is embedded in an environment
#' within the RTMB object.
#'
#' @param MSAassess [MSAassess-class] object returned by `fit_MSA()`
#' @return [MSAdata-class] object
#' @export
get_MSAdata <- function(MSAassess) {
  func <- attr(MSAassess@obj$env$data, "func")
  version <- strsplit(attr(MSAassess, "version"), "multiSA ")[[1]][2]

  if (version >= "0.2.0") {
    d <- "x"
  } else {
    d <- "MSAdata"
  }
  MSAdata <- get(d, envir = environment(func), inherits = FALSE)
  return(MSAdata)
}

make_unique_names <- function(x, select = c("fixed", "random", "report")) {
  select <- match.arg(select)

  if (select == "fixed") {
    par_names <- unique(names(x$par.fixed))
    par_list <- as.list(x, what = "Estimate", report = FALSE)
  } else if (select == "report") {
    par_names <- unique(names(x$value))
    par_list <- as.list(x, what = "Estimate", report = TRUE)
  } else {
    par_names <- unique(names(x$random))
    par_list <- as.list(x, what = "Estimate", report = FALSE)
  }


  par_dims <- lapply(par_names, function(y) {
    dim_y <- dim(par_list[[y]])
    if (is.null(dim_y)) dim_y <- length(par_list[[y]])
    ind_y <- lapply(dim_y, function(i) seq(1, i))

    est_grid <- do.call(expand.grid, ind_y)

    if (select != "report") {
      map_y <- attr(x$env$parameters[[y]], "map")

      if (!is.null(map_y)) {
        est_y <- map_y >= 0 & !duplicated(map_y)
        est_grid <- est_grid[est_y, ]
      }
    }

    dim_char <- sapply(1:nrow(est_grid), function(i) {
      paste0(
        "[",
        paste0(est_grid[i, ], collapse = ", "),
        "]"
      )
    })
    paste0(y, dim_char)
  })

  do.call(c, par_dims)
}


make_yearseason <- function(year, nm = 4) {
  if (nm <= 1) return(year)
  year_long <- lapply(year, function(y) y + (1:nm - 1)/nm)
  do.call(c, year_long)
}

#' @importFrom reshape2 acast melt
collapse_yearseason <- function(x) {
  dim_x <- dim(x)
  if (length(dim_x) > 2) {
    dimnames(x)[[1]] <- 1:dim_x[1]
    dimnames(x)[[2]] <- 1:dim_x[2]

    x_df <- reshape2::melt(x)
    x_df$Y <- as.numeric(x_df$Var1) + (as.numeric(x_df$Var2) - 1)/dim_x[2]

    dims <- c("Y", paste0("Var", 3:length(dim_x))) |> as.list()
    xout <- reshape2::acast(x_df, dims, value.var = "value")
    dimnames(xout) <- NULL

    return(xout)

  } else {
    return(as.numeric(t(x)))
  }
}

message <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) |> paste0(collapse = "") |> usethis::ui_done()
  } else {
    base::message(...)
  }
}


message_info <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) |> paste0(collapse = "") |> usethis::ui_info()
  } else {
    base::message(...)
  }
}

message_oops <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) |> paste0(collapse = "") |> usethis::ui_oops()
  } else {
    base::message(...)
  }
}

warning <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) |> paste0(collapse = "") |> usethis::ui_warn()
  } else {
    base::warning(...)
  }
}


stop <- function(..., call. = TRUE, domain = NULL) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) |> paste0(collapse = "") |> usethis::ui_stop()
  } else {
    base::stop(..., call. = call., domain = domain)
  }
}

