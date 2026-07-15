
#' Selectivity at age and length
#'
#' @description
#' Calculate selectivity at age and length from a matrix of parameters.
#' * [conv_selpar()] converts parameters from log or logit space to real units.
#' * [calc_sel_len()] calculates selectivity at length.
#' * [calc_fsel_age()] calculates selectivity at age for fisheries, and coordinates dummy fleets.
#' * [calc_isel_age()] calculates selectivity at age for indices, and can map selectivity from fisheries
#' or population parameters (e.g, mature or total biomass).
#'
#' @param x Estimated parameters. Matrix `[3, f]`
#' @param type,fsel_type Character string to indicate the functional form of selectivity. See details below.
#' @param maxage Maximum value of the age of full selectivity
#' @param maxL Maximum value of the length of full selectivity
#' @details
#' Options for argument `type` include:
#'
#' - functional forms with respect to length: `"logistic_length", "dome_length"`
#' - functional forms with respect to age: `"logistic_age", "dome_age"`
#' - for surveys, an integer (`f`) to map index selectivity at age to fleet `f` (will be coerced to integer)
#' - `"mature"` to fix to maturity at age schedule (`"SB"` is also acceptable for backwards compatibility)
#' - `"total"` to fix selectivity to 1 for all ages (`"B"` is also acceptable for backwards compatibility)
#' - `"length_x_y"` to specify selectivity to 1 between lengths `x` and `y` (example: `"length_20_50"`)
#' - `"age_x_y"` to specify selectivity to 1 between age `x` and `y` (example: `"age_2_4"`)
#' - for surveys, `"f_x_y"` uses selectivity values from fleet `f` between bins `x` and `y`
#' (either length or age depending on definition of selectivity for `f`) (example: `"3_2_4"`)
#'
#' @section Converting selectivity parameters (conv_selpar):
#' The first row of `x` corresponds to the length or age of full selectivity: \eqn{\mu_f = L_{max}/(1 + \exp(-x_{1,f}))}
#'
#' The second row of `x` corresponds to the width of the ascending limb for selectivity:
#' \eqn{\sigma_f^{asc} = \exp(x_{2,f})}
#'
#' The third row of `x` corresponds to the width of the descending limb for selectivity (if dome-shaped):
#' \eqn{\sigma_f^{des} = \exp(x_{3,f})}
#' @return
#' [conv_selpar()] returns a matrix of the same dimensions as `x`.
#' @export
conv_selpar <- function(x, type, maxage, maxL) {
  nf <- length(type)
  stopifnot(ncol(x) == nf)
  stopifnot(nrow(x) >= 3)

  sel_par <- sapply(1:nf, function(f) {
    sd_asc <- exp(x[2, f])
    sd_desc <- exp(x[3, f])
    if (grepl("logistic_age|dome_age", type[f])) {
      Aapical <- maxage * plogis(x[1, f])
      v <- c(Aapical, sd_asc, sd_desc)
    } else if (grepl("logistic_length|dome_length", type[f])) {
      Lapical <- maxL * plogis(x[1, f])
      v <- c(Lapical, sd_asc, sd_desc)
    } else {
      v <- rep(NA_real_, 3)
    }
    return(v)
  })

  return(sel_par)
}

#' @rdname conv_selpar
#' @aliases calc_sel_len
#' @param sel_par Matrix of parameters returned by [conv_selpar()]
#' @param lmid Midpoints of length bins for calculating selectivity at length
#' @param fsel_len Selectivity at length matrix for fleets, returned by previous call to `calc_sel_len()`
#' @section Length selectivity (calc_sel_len):
#' The selectivity at length is
#' \deqn{
#' s_{\ell} =
#' \begin{cases}
#' \exp(\alpha) & L_{\ell} < \mu_f\\
#' \exp(\beta) & L_{\ell} \ge \mu_f\\
#' \end{cases}
#' }
#' where
#' \eqn{
#' \alpha = -0.5(L_\ell - \mu_f)^2/(\sigma_f^{asc})^2
#' }
#' and
#' \eqn{
#' \beta = -0.5(L_\ell - \mu_f)^2/(\sigma_f^{des})^2
#' }
#' @return
#' [calc_sel_len()] returns a matrix `[l, f]`, i.e., `[length(lmid), length(type)]`.
#' @export
calc_sel_len <- function(sel_par, lmid, type, fsel_type, fsel_len) {
  nf <- length(type)
  is_ad <- inherits(sel_par, "advector")

  sel_lf <- sapply(1:nf, function(f) {
    if (grepl("logistic_length|dome_length", type[f])) {
      ex_asc <- (lmid - sel_par[1, f])/sel_par[2, f]
      asc <- exp(-0.5 * ex_asc * ex_asc)

      if (grepl("logistic_length", type[f])) {
        desc <- 1
      } else {
        ex_desc <- (lmid - sel_par[1, f])/sel_par[3, f]
        desc <- exp(-0.5 * ex_desc * ex_desc)
      }
      v <- CondExpLt(lmid, sel_par[1, f], asc, desc)
      v <- v/max(v)
    } else if (startsWith(type[f], "length")) {
      sel_char <- strsplit(type[f], "_")[[1]]
      lmin <- as.numeric(sel_char[2])
      lmax <- as.numeric(sel_char[3])
      v <- ifelse(lmid >= lmin & lmid <= lmax, 1, 0)
    } else {
      sel_char <- strsplit(type[f], "_")[[1]]

      if (length(sel_char) == 3) {
        ff <- suppressWarnings(as.numeric(sel_char)[1])

        if (is.numeric(ff) && grepl("length", fsel_type[ff])) {
          lmin <- as.numeric(sel_char[2])
          lmax <- as.numeric(sel_char[3])
          v <- ifelse(lmid >= lmin & lmid <= lmax, fsel_len[, ff], 0)
        } else {
          v <- rep(NA_real_, length(lmid))
        }
      } else {
        v <- rep(NA_real_, length(lmid))
      }
    }
    if (is_ad) v <- advector(v)
    return(v)
  })

  return(sel_lf)
}

#' @rdname conv_selpar
#' @param sel_len Selectivity at length matrix returned by [calc_sel_len()]
#' @param LAK Length-at-age probability matrix. Matrix `[a, length(lmid)]`
#' @param sel_block Integer vector. Length `length(type)`. See details below.
#' @param mat Maturity at age vector
#' @param a Integer vector of ages corresponding to the rows of `LAK` (as well as `mat`)
#' @section Age selectivity (calc_fsel_age):
#' The equivalent selectivity at age is converted from the length values (`sel_len`) as
#' \deqn{
#' s_a = \sum_\ell s_\ell \times \textrm{Prob}(L_{\ell}|a)
#' }
#'
#' If selectivity is explicitly in age units, then the values are directly calculated
#' from parameters in `sel_par`.
#'
#' Vector `sel_block` assigns the output selectivity from a different column of the input matrix
#' and facilitates time-varying selectivity in time blocks. For example, `sel_block[1] <- 2` means
#' that selectivity values in the first column of the output is based on the second column of the
#' input matrices (`sel_len[, 2]` or `sel_par[, 2]`).
#' @return
#' [calc_fsel_age()] returns a matrix `[a, f]`, i.e., `[a, length(sel_block)]`
#' @export
calc_fsel_age <- function(sel_len, LAK, type, sel_par, sel_block = seq(1, length(type)), mat, a = seq(1, nrow(LAK)) - 1) {
  nf <- length(sel_block)

  is_ad <- inherits(sel_par, "advector")

  sel_af <- sapply(1:nf, function(ff) {
    f <- sel_block[ff]
    if (grepl("length", type[f])) {
      v <- sel_len[, f] %*% t(LAK)
      #v <- v/max(v)
    } else if (grepl("logistic_age|dome_age", type[f])) {

      ex_asc <- (a - sel_par[1, f])/sel_par[2, f]
      asc <- exp(-0.5 * ex_asc * ex_asc)

      if (grepl("logistic_age", type[f])) {
        desc <- 1
      } else {
        ex_desc <- (a - sel_par[1, f])/sel_par[3, f]
        desc <- exp(-0.5 * ex_desc * ex_desc)
      }
      v <- CondExpLt(a, sel_par[1, f], asc, desc)
      v <- v/max(v)
    } else if (startsWith(type[f], "age")) {
      sel_char <- strsplit(type[f], "_")[[1]]
      amin <- as.numeric(sel_char[2])
      amax <- as.numeric(sel_char[3])
      v <- ifelse(a >= amin & a <= amax, 1, 0)
    } else if (type[f] %in% c("SB", "mature")) {
      v <- mat
    } else if (type[f] %in% c("B", "total")) {
      v <- rep(1, length(a))
    } else if (type[f] == "free") {
      v <- plogis(sel_par[, f])
    }
    if (is_ad) v <- advector(v)
    return(v)
  })

  return(sel_af)
}

#' @rdname conv_selpar
#' @param fsel_age Matrix returned by [calc_fsel_age()]
#' @return
#' [calc_isel_age()] returns a matrix `[a, i]`, i.e., `[a, length(type)]`
#' @export
calc_isel_age <- function(sel_len, LAK, type, sel_par, fsel_age, maxage, mat,
                          a = seq(1, nrow(LAK)) - 1,
                          fsel_type, fsel_len) {
  ni <- length(type)
  is_ad <- inherits(sel_par, "advector") || inherits(fsel_age, "advector")

  sel_ai <- sapply(1:ni, function(i) {
    ti <- type[i]
    tii <- suppressWarnings(as.integer(ti))

    if (is.na(tii)) { # Not an integer
      if (grepl("length", ti)) {
        v <- sel_len[, i] %*% t(LAK)
      } else if (grepl("age", ti)) {

        if (grepl("logistic_age|dome_age", ti)) {
          ex_asc <- (a - sel_par[1, i])/sel_par[2, i]
          asc <- exp(-0.5 * ex_asc * ex_asc)

          if (grepl("logistic_age", ti)) {
            desc <- 1
          } else {
            ex_desc <- (a - sel_par[1, i])/sel_par[3, i]
            desc <- exp(-0.5 * ex_desc * ex_desc)
          }
          v <- CondExpLt(a, sel_par[1, i], asc, desc)
          v <- v/max(v)

        } else if (startsWith(ti, "age")) {

          sel_char <- strsplit(ti, "_")[[1]]
          amin <- as.numeric(sel_char[2])
          amax <- as.numeric(sel_char[3])
          v <- ifelse(a >= amin & a <= amax, 1, 0)

        }

      } else if (ti == "SB") {
        v <- mat
      } else if (ti == "B") {
        v <- rep(1, length(a))
      } else if (ti == "free") {
        v <- plogis(sel_par[, i])
      } else {
        sel_char <- strsplit(ti, "_")[[1]]
        ff <- suppressWarnings(as.numeric(sel_char[1]))
        sel_f <- fsel_type[ff]
        if (grepl("length", sel_f)) {
          v <- sel_len[, i] %*% t(LAK)
        } else {
          stop("Error: can't figure out how to calculate index selectivity")
        }
      }
    } else if (is.integer(tii)) {
      v <- fsel_age[, tii]
    } else {
      stop("Error: can't identify fleet number to mirror index selectivity")
    }
    if (is_ad) v <- advector(v)
    return(v)
  })

  return(sel_ai)
}
