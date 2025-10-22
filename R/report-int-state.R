

#' @importFrom grDevices hcl.colors
make_color <- function(n, type = c("fleet", "region", "stock"), alpha = 1) {
  if (n > 1) {
    type <- match.arg(type)
    pal <- switch(
      type,
      "fleet" = "Set2",
      "region" = "Sunset",
      "stock" = "Plasma"
    )
    grDevices::hcl.colors(n, palette = pal, alpha = alpha)
  } else {
    "black"
  }
}

#' @importFrom tinyplot tinyplot
make_tinyplot <- function(year, x, ylab, name, color, type = "o",
                          leg = if (ncol(x) > 1) substitute(legend(title = NULL)) else "none",
                          ylim = c(0, 1.1) * range(x, na.rm = TRUE)) {

  df <- structure(x, dimnames = list(year = year, by = name)) %>%
    reshape2::melt()

  tinyplot_args <- list(
    x = df$year, y = df$value, by = factor(df$by), xlab = "Year", ylab = ylab,
    type = type,
    legend = leg,
    col = color,
    pch = 16,
    ylim = ylim,
    grid = TRUE
  )
  do.call(tinyplot, tinyplot_args)
  abline(h = 0, col = "grey60")
  box()

  invisible()
}

#' @importFrom graphics box legend
#' @importFrom tinyplot tinyplot type_barplot type_lines
#' @importFrom reshape2 melt
barplot2 <- function(x, cols, leg.names, facet.names = NULL, xval, ylab = ifelse(prop, "Proportion", "Value"),
                     border = ifelse(nrow(x) > 60, NA, "grey60"), prop = TRUE, facet.free = FALSE) {

  if (length(x) == dim(x)[1]) {
    if (prop) x[] <- 1
    plot(xval, as.numeric(x), xlab = "Year", ylab = ylab, ylim = c(0, 1.1) * range(x), pch = 16, type = "o", zero_line = TRUE)
    return(invisible())
  }

  ndim <- length(dim(x))
  if (ndim == 2) {
    if (prop) {
      p <- apply(x, 1, function(x) x/sum(x, na.rm = TRUE))
      if (is.null(dim(p))) p <- matrix(p, 1, length(xval))
      ylim <- c(0, 1)
    } else {
      p <- t(x)
      ylim <- c(0, 1.1) * range(rowSums(x))
    }

    df <- structure(p, dimnames = list(by = leg.names, xval = xval)) %>%
      reshape2::melt()

  } else if (ndim == 3) {
    if (prop) {
      p <- apply(x, c(1, 3), function(x) x/sum(x, na.rm = TRUE))
      if (is.null(dim(p))) p <- matrix(p, 1, length(xval))
      ylim <- c(0, 1)
    } else {
      p <- aperm(x, c(2, 1, 3))
      ylim <- c(0, 1.1) * range(apply(p, 1:2, sum))
    }

    df <- structure(p, dimnames = list(by = leg.names, xval = xval, facet = facet.names)) %>%
      reshape2::melt()

  } else {
    stop("Dimension not usable by barplot2")
  }

  do.facet <- ndim == 3 && length(unique(df$facet)) > 1
  do.group <- dim(p)[1] > 1

  if (missing(cols)) {
    ncat <- nrow(p)
    cols <- make_color(ncat)
  }

  if (do.group) {

    tinyplot_args <- list(
      x = df$xval, y = df$value,
      by = df$by,
      grid = TRUE,
      legend = substitute(legend(title = NULL)),
      facet = if (do.facet) df$facet else NULL,
      facet.args = if (do.facet) list(free = facet.free) else NULL,
      xlab = "Year", ylab = ylab,
      type = type_barplot(width = 1),
      palette = cols,
      col = border,
      ylim = NULL
    )

  } else {

    tinyplot_args <- list(
      x = df$xval, y = df$value,
      by = NULL,
      grid = TRUE,
      legend = NULL,
      facet = if (do.facet) df$facet else NULL,
      facet.args = if (do.facet) list(free = facet.free) else NULL,
      xlab = "Year", ylab = ylab,
      type = type_lines(),
      palette = NULL,
      col = cols,
      ylim = c(0, 1.1) * range(df$value)
    )

  }

  do.call(tinyplot, tinyplot_args)
  if (!do.facet) box()

  invisible()
}


# State variable plots ----
#' @name plot-MSA-state
#' @title Plotting functions for fitted MSA model
#' @description A set of functions to plot state variables (biomass, recruitment time series, etc.)
#' @return Various base graphics plots
NULL

#' @rdname plot-MSA-state
#' @aliases plot_S
#'
#' @param fit [MSAassess-class] object returned by [fit_MSA()]
#' @param by Character to indicate dimension for multivariate plots
#' @param s Integer for the corresponding stock
#' @param prop Logical, whether to plot proportions (TRUE) or absolute numbers
#' @param facet_free Logical, whether to allow the y-axis limits to vary by panel in facetted plots
#' @details
#' - `plot_S` plots spawning output by stock or region (either whole numbers or proportions for the latter)
#'
#' @export
plot_S <- function(fit, by = c("stock", "region"), r, s, prop = FALSE, facet_free = FALSE) {
  by <- match.arg(by)
  var <- "S_yrs"

  d <- get_MSAdata(fit)
  Dlabel <- d@Dlabel
  Dmodel <- d@Dmodel

  year <- Dlabel@year
  ny <- length(year)

  if (missing(r)) r <- 1:Dmodel@nr
  if (missing(s)) s <- 1:Dmodel@ns
  rname <- Dlabel@region[r]
  sname <- Dlabel@stock[s]

  if (by == "stock") {
    leg.name <- sname
    facet.name <- rname
    x <- array(fit@report[[var]][, r, s, drop = FALSE], c(ny, length(rname), length(sname))) %>%
      aperm(c(1, 3, 2))
  } else {
    leg.name <- rname
    facet.name <- sname
    x <- array(fit@report[[var]][, r, s, drop = FALSE], c(ny, length(rname), length(sname)))
  }

  color <- make_color(ncol(x), type = by)

  if (prop) {
    ylab <- "Spawning fraction"
  } else {
    ylab <- "Spawning output"
  }

  barplot2(x, cols = color, leg.names = leg.name, facet.names = facet.name, xval = year, ylab = ylab, prop = prop,
           facet.free = facet_free)

  invisible(array2DF(x, responseName = "S"))
}


#' @rdname plot-MSA-state
#' @aliases plot_B
#' @details
#' - `plot_B` plots total biomass by stock or region (either whole numbers or proportions for the latter)
#'
#' @export
plot_B <- function(fit, by = c("stock", "region"), r, s, prop = FALSE, facet_free = FALSE) {
  by <- match.arg(by)
  var <- "B_ymrs"

  d <- get_MSAdata(fit)
  Dlabel <- d@Dlabel
  Dmodel <- d@Dmodel

  year <- Dlabel@year
  ny <- length(year)
  nm <- max(length(Dlabel@season), 1)

  if (missing(r)) r <- 1:Dmodel@nr
  if (missing(s)) s <- 1:Dmodel@ns
  rname <- Dlabel@region[r]
  sname <- Dlabel@stock[s]

  if (by == "stock") {
    leg.name <- sname
    facet.name <- rname
    x <- array(fit@report[[var]][, , r, s, drop = FALSE], c(ny, nm, length(rname), length(sname))) %>%
      aperm(c(1, 2, 4, 3)) # B_ymsr
  } else {
    leg.name <- rname
    facet.name <- sname
    x <- array(fit@report[[var]][, , r, s, drop = FALSE], c(ny, nm, length(rname), length(sname))) # B_ymrs
  }

  year <- make_yearseason(year, nm)
  x <- collapse_yearseason(x)

  color <- make_color(ncol(x), type = by)

  if (prop) {
    ylab <- "Biomass fraction"
  } else {
    ylab <- "Total biomass"
  }

  barplot2(x, cols = color, leg.names = leg.name, facet.names = facet.name, xval = year, ylab = ylab, prop = prop,
           facet.free = facet_free)

  invisible(array2DF(x, responseName = "B"))
}



#' @rdname plot-MSA-state
#' @aliases plot_R
#' @details
#' - `plot_R` plots recruitment by stock
#'
#' @export
plot_R <- function(fit, s) {
  var <- "R_ys"

  Dlabel <- get_MSAdata(fit)@Dlabel
  year <- Dlabel@year

  if (missing(s)) {
    x <- fit@report[[var]]
    name <- Dlabel@stock
    ylab <- "Recruitment"
  } else {
    x <- fit@report[[var]][, s, drop = FALSE]
    name <- Dlabel@stock[s]
    ylab <- paste(name, "recruitment")
  }

  color <- make_color(ncol(x), "stock")
  make_tinyplot(year, x, ylab, name, color)
  #matplot(year, x, xlab = "Year", ylab = ylab, type = "o", col = color, pch = 16, lty = 1,
  #        ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
  #if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  x <- structure(x, dimnames = list(year = year, stock = name))
  invisible(array2DF(x, "R"))
}


#' @rdname plot-MSA-state
#' @aliases plot_SRR
#' @param phi Logical, whether to plot unfished replacement line
#' @details
#' - `plot_SRR` plots the stock-recruitment relationship and history (realized recruitment) by stock
#' @importFrom graphics points lines
#' @export
plot_SRR <- function(fit, s = 1, phi = TRUE) {
  dat <- get_MSAdata(fit)
  Dlabel <- get_MSAdata(fit)@Dlabel

  S_y <- apply(fit@report$S_yrs[, , s, drop = FALSE], 1, sum)
  R_y <- fit@report$R_ys[, s]

  Rpred_y <- R_y/fit@report$Rdev_ys[, s]

  a <- fit@report$sralpha_s[s]
  b <- fit@report$srbeta_s[s]

  S_SRR <- seq(0, max(S_y), length.out = 100)
  R_SRR <- calc_recruitment(S_SRR, SRR = dat@Dstock@SRR_s[s], a = a, b = b)

  S2 <- S_y[-1]
  R2 <- Rpred_y[-1]

  plot(S_SRR, R_SRR, type = "l", lwd = 1, lty = 3,
       xlab = "Spawning output", ylab = "Recruitment",
       xaxs = "i", yaxs = "i", xlim = c(0, 1.1) * range(S_y), ylim = c(0, 1.1) * range(R_y))
  lines(S2[order(S2)], R2[order(R2)], lwd = 2)
  points(S_y, R_y, pch = 16)

  if (phi) {
    phi_s <- fit@report$phi_s[s]
    abline(a = 0, b = 1/phi_s, lty = 2)
  }
  invisible()
}


#' @rdname plot-MSA-state
#' @aliases plot_Rdev
#' @param log Logical, whether to plot the natural logarithm of the response variable
#' @details
#' - `plot_Rdev` plots recruitment deviations by stock
#' @importFrom graphics arrows
#' @export
plot_Rdev <- function(fit, s = 1, log = TRUE) {

  Dlabel <- get_MSAdata(fit)@Dlabel
  year <- Dlabel@year

  if (log) {
    if (length(fit@SD) > 1) {
      x <- as.list(fit@SD, what = "Estimate")$log_rdev_ys[, s]
      std <- as.list(fit@SD, what = "Std. Error")$log_rdev_ys[, s]
      std[is.na(std)] <- 0
    } else {
      x <- fit@obj$env$parList(par = fit@obj$env$last.par.best)$log_rdev_ys[, s]
      std <- numeric(length(x))
    }

    upper <- x + 1.96 * std
    lower <- x - 1.96 * std

    plot(year, x, xlab = "Year", ylab = "log Recruitment deviations", type = "o", pch = 16,
         ylim = range(lower, upper), lty = 3)
    arrows(x0 = year, y0 = lower, y1 = upper, length = 0)
    abline(h = 0, lty = 2)

  } else {
    x <- fit@report$Rdev_ys[, s]
    plot(year, x, xlab = "Year", ylab = "Recruitment deviations", type = "o", pch = 16,
         ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)

    abline(h = 1, lty = 2)
  }
  invisible()
}

#' @rdname plot-MSA-state
#' @aliases plot_Fstock
#' @details
#' - `plot_Fstock` plots apical instantaneous fishing mortality (per year or per season) by stock
#'
#' @export
plot_Fstock <- function(fit, s, by = c("annual", "season")) {
  by <- match.arg(by)

  dat <- get_MSAdata(fit)
  year <- dat@Dlabel@year
  nm <- dat@Dmodel@nm

  unit <- ifelse(by == "annual", "year", "season")

  if (by == "annual") {
    var <- "F_yas"
    if (missing(s)) {
      x <- fit@report[[var]]
    } else {
      x <- fit@report[[var]][, , s, drop = FALSE]
    }
    x <- apply(x, c(1, 3), max)

  } else if (by == "season" && nm > 1) {
    if (missing(s)) {
      F_ymas <- sapply2(1:dat@Dmodel@ns, function(s) {
        sapply2(1:dat@Dmodel@na, function(a) {
          sapply(1:nm, function(m) {
            sapply(1:dat@Dmodel@ny, function(y) {
              N <- sum(fit@report$N_ymars[y, m, a, , s])
              CN <- sum(fit@report$CN_ymafrs[y, m, a, , , s])
              calc_summary_F(M = fit@report$M_yas[y, a, s]/nm, N = N, CN = CN, Fmax = 100)
            })
          })
        })
      })
    } else {
      F_ymas <- sapply2(1, function(...) {
        sapply2(1:dat@Dmodel@na, function(a) {
          sapply(1:nm, function(m) {
            sapply(1:dat@Dmodel@ny, function(y) {
              N <- sum(fit@report$N_ymars[y, m, a, , s])
              CN <- sum(fit@report$CN_ymafrs[y, m, a, , , s])
              calc_summary_F(M = fit@report$M_yas[y, a, s]/nm, N = N, CN = CN, Fmax = 100)
            })
          })
        })
      })
    }
    x <- apply(F_ymas, c(1, 2, 4), max) %>%
      collapse_yearseason()

    year <- make_yearseason(year, nm)
  }

  if (exists("x", inherits = FALSE)) {
    x[is.infinite(x)] <- NA

    if (missing(s)) {
      name <- dat@Dlabel@stock
      ylab <- paste0("Apical fishing mortality (per ", unit, ")")
    } else {
      name <- dat@Dlabel@stock[s]
      ylab <- paste0(name, " apical fishing mortality (per ", unit, ")")
    }

    color <- make_color(ncol(x), "stock")
    make_tinyplot(year, x, ylab, name, color)
    #matplot(year, x, xlab = "Year", ylab = ylab, type = "o", col = color, pch = 16, lty = 1,
    #        ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
    #if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  }

  invisible()
}



#' @rdname plot-MSA-state
#' @aliases plot_self
#' @param f Integer for the corresponding fleet
#' @param type For `plot_self`, indicates whether to plot the selectivity by age or length.
#' @details
#' - `plot_self` plots fishery selectivity
#' @export
plot_self <- function(fit, f = 1, type = c("length", "age")) {
  type <- match.arg(type)

  dat <- get_MSAdata(fit)
  Dfishery <- dat@Dfishery

  sel_block <- Dfishery@sel_block_yf[, f]
  sel_b <- Dfishery@sel_f[unique(sel_block)]
  fname <- dat@Dlabel@fleet[f]

  year <- dat@Dlabel@year

  if (type == "length" && all(grepl("length", sel_b))) {
    lmid <- dat@Dmodel@lmid
    x <- fit@report$sel_lf[, unique(sel_block), drop = FALSE]

    color <- make_color(ncol(x), "fleet")
    matplot(lmid, x, xlab = "Length", ylab = paste(fname, "selectivity"),
            type = "o", col = color, pch = 16,
            ylim = c(0, 1), lty = 1, zero_line = TRUE)
    if (ncol(x) > 1) {
      name <- sapply(unique(sel_block), function(i) {
        y <- year[sel_block == i]
        if (length(y) == 1) {
          return(y)
        } else {
          return(paste(range(y), collapse = "-"))
        }
      })
      legend("topright", legend = name, col = color, lwd = 1, pch = 16)
    }
  } else if (type == "age") {

    m <- 1
    s <- 1

    x <- fit@report$sel_ymafs[, m, , f, s]
    xx <- apply(x, 2, diff)

    if (any(xx != 0)) {
      ybreak <- c(1, which(rowSums(xx) > 0) + 1)
      name <- sapply(1:length(ybreak), function(i) {
        if (i == length(ybreak)) {
          y <- c(year[ybreak[i]], year[length(year)])
        } else {
          y <- year[c(ybreak[i], ybreak[i+1] - 1)]
        }
        paste(range(y), collapse = "-")
      })
      x <- x[ybreak, , drop = FALSE]

    } else {
      x <- x[1, , drop = FALSE]
    }
    age <- dat@Dlabel@age

    color <- make_color(nrow(x), "fleet")
    matplot(age, t(x), xlab = "Age", ylab = paste(fname, "selectivity"),
            type = "o", col = color, pch = 16,
            ylim = c(0, 1), lty = 1, zero_line = TRUE)
    if (nrow(x) > 1) legend("topright", legend = name, col = color, lwd = 1, pch = 16)
  }
  invisible()
}

#' @rdname plot-MSA-state
#' @aliases plot_seli
#' @param i Integer for the corresponding survey
#' @details
#' - `plot_seli` plots index selectivity
#' @export
plot_seli <- function(fit, i = 1) {
  dat <- get_MSAdata(fit)
  sel_i <- dat@Dsurvey@sel_i[i]
  mirror_f <- suppressWarnings(as.numeric(sel_i))

  iname <- dat@Dlabel@index[i]

  if (!is.na(mirror_f)) {
    plot_self(fit, f = mirror_f)
  } else if (grepl("length", sel_i)) {
    lmid <- dat@Dmodel@lmid
    x <- fit@report$sel_li[, i]

    plot(lmid, x, xlab = "Length", ylab = paste(iname, "selectivity"),
         type = "o", pch = 16,
         ylim = c(0, 1), lty = 1, zero_line = TRUE)
  } else {

    m <- 1
    s <- 1

    x <- fit@report$sel_ymais[, m, , i, s]
    xx <- apply(x, 2, diff)

    if (any(xx != 0)) {
      year <- dat@Dlabel@year
      ybreak <- c(1, which(rowSums(xx) > 0) + 1)
      name <- sapply(1:length(ybreak), function(i) {
        if (i == length(ybreak)) {
          y <- c(year[ybreak[i]], year[length(year)])
        } else {
          y <- year[c(ybreak[i], ybreak[i+1] - 1)]
        }
        paste(range(y), collapse = "-")
      })
      x <- x[ybreak, , drop = FALSE]

    } else {
      x <- x[1, , drop = FALSE]
    }
    age <- dat@Dlabel@age
    color <- make_color(nrow(x), "fleet")

    matplot(age, t(x), xlab = "Age", ylab = paste(iname, "selectivity"),
            type = "o", col = color, pch = 16,
            ylim = c(0, 1), lty = 1, zero_line = TRUE)
    if (nrow(x) > 1) legend("topright", legend = name, col = color, lwd = 1, pch = 16)
  }
  invisible()
}


#' @rdname plot-MSA-state
#' @aliases plot_selstock
#' @param plot2d Character, plotting function for either a [contour()] or [filled.contour()] plot
#' @param by Character to indicate whether to calculate selectivity from F per year or per season
#' @param ... Other arguments to the base graphics function
#' @details
#' - `plot_selstock` plots the realized selectivity from total catch and total abundance at age
#' @export
#' @importFrom graphics contour filled.contour
plot_selstock <- function(fit, s = 1, by = c("annual", "season"), plot2d = c("contour", "filled.contour"), ...) {

  by <- match.arg(by)

  plot2d <- match.arg(plot2d)
  plot2d <- match.fun(plot2d)

  dat <- get_MSAdata(fit)
  year <- dat@Dlabel@year
  age <- dat@Dlabel@age

  nm <- max(length(dat@Dlabel@season), 1)

  if (by == "annual") {
    sel_ya <- fit@report$F_yas[, , s] %>%
      apply(1, function(x) x/max(x)) %>%
      t()
  } else if (by == "season" && nm > 1) {

    year <- make_yearseason(year, nm)
    F_yma <- sapply2(1:dat@Dmodel@na, function(a) {
      sapply(1:nm, function(m) {
        sapply(1:dat@Dmodel@ny, function(y) {
          N <- sum(fit@report$N_ymars[y, m, a, , s])
          CN <- sum(fit@report$CN_ymafrs[y, m, a, , , s])
          calc_summary_F(M = fit@report$M_yas[y, a, s]/nm, N = N, CN = CN, Fmax = 100)
        })
      })
    })
    sel_ya <- collapse_yearseason(F_yma) %>%
      apply(1, function(x) x/max(x)) %>%
      t()
  }

  if (exists("sel_ya", inherits = FALSE)) {
    plot2d(x = year, y = age, xlab = "Year", ylab = "Age", z = sel_ya, levels = seq(0, 1, 0.1), ...)
  }

  invisible()
}

#' @rdname plot-MSA-state
#' @aliases plot_N
#' @param m Integer for the corresponding season
#' @param r Integer for the corresponding region
#' @param ... Other argument to the base graphics function
#' @details
#' - `plot_N` reports total abundance at age
#' @export
#' @importFrom graphics contour filled.contour
plot_N <- function(fit, m = 1, r, s = 1, plot2d = c("contour", "filled.contour"), ...) {
  plot2d <- match.arg(plot2d)
  plot2d <- match.fun(plot2d)

  dat <- get_MSAdata(fit)
  if (missing(r)) r <- 1:dat@Dmodel@nr
  if (length(m) > 1) stop("length(m) should be one")
  ny <- dat@Dmodel@ny
  year <- dat@Dlabel@year
  age <- dat@Dlabel@age

  N_ya <- apply(fit@report$N_ymars[1:ny, m, , r, s, drop = FALSE], c(1, 3), sum)

  dots <- list(...)
  if (!length(dots$nlevels)) dots$nlevels <- 10
  if (!length(dots$levels)) dots$levels <- exp(pretty(log(range(N_ya)), 10))

  dots$x <- year
  dots$y <- age
  dots$xlab = "Year"
  dots$ylab <- "Age"
  dots$z <- N_ya

  do.call(plot2d, dots)

  invisible()
}


#' @rdname plot-MSA-state
#' @aliases plot_V
#' @details
#' - `plot_V` plots vulnerable biomass, availability to the fishery
#' @export
plot_V <- function(fit, f = 1, by = c("stock", "region"), prop = FALSE, facet_free = FALSE) {
  by <- match.arg(by)
  var <- "VB_ymfrs"

  d <- get_MSAdata(fit)
  Dlabel <- d@Dlabel
  Dmodel <- d@Dmodel

  year <- Dlabel@year
  ny <- length(year)
  nm <- max(length(Dlabel@season), 1)

  r <- 1:Dmodel@nr
  s <- 1:Dmodel@ns
  rname <- Dlabel@region[r]
  sname <- Dlabel@stock[s]

  if (by == "stock") {
    leg.name <- sname
    facet.name <- rname
    x <- array(fit@report[[var]][, , f, , , drop = FALSE], c(ny, nm, length(rname), length(sname))) %>%
      aperm(c(1, 2, 4, 3)) # B_ymsr
  } else {
    leg.name <- rname
    facet.name <- sname
    x <- array(fit@report[[var]][, , f, , , drop = FALSE], c(ny, nm, length(rname), length(sname))) # B_ymrs
  }

  year <- make_yearseason(year, nm)
  x <- collapse_yearseason(x)

  color <- make_color(ncol(x), type = by)

  fname <- Dlabel@fleet[f]

  if (prop) {
    ylab <- paste("Proportion biomass available to", fname)
  } else {
    ylab <- paste("Biomass available to", fname)
  }

  barplot2(x, cols = color, leg.names = leg.name, facet.names = facet.name, xval = year, ylab = ylab, prop = prop,
           facet.free = facet_free)

  invisible(array2DF(x, responseName = "V"))
}


#' @rdname plot-MSA-state
#' @aliases plot_Ffleet
#' @details
#' - `plot_Ffleet` plots apical instantaneous fishing mortality (per season) by fleet
#'
#' @export
plot_Ffleet <- function(fit, f = 1) {
  var <- "F_ymfr"

  Dlabel <- get_MSAdata(fit)@Dlabel
  year <- Dlabel@year
  nm <- max(length(Dlabel@season), 1)

  x <- apply(fit@report[[var]][, , f, , drop = FALSE], c(1, 2, 4), identity)
  name <- Dlabel@region

  year <- make_yearseason(year, nm)
  x <- collapse_yearseason(x)

  color <- make_color(ncol(x), type = "region")

  fname <- Dlabel@fleet[f]
  ylab <- paste(fname, "fishing mortality")

  make_tinyplot(year, x, ylab, name, color)
  #matplot(year, x, xlab = "Year", ylab = ylab, type = "o", col = color, pch = 16, lty = 1,
  #        ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
  #if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  invisible()





}

#' @rdname plot-MSA-state
#' @aliases plot_mov
#' @param y Integer, year for plotting the movement matrix (last model year is the default)
#' @param a Integer, corresponding age for plotting the movement matrix (age 1 is the default)
#' @param palette Character, palette name to pass to [grDevices::hcl.colors()]. See [grDevices::hcl.pals()] for options.
#' @details
#' - `plot_mov` plots movement matrices and the corresponding equilibrium distribution in multi-area models
#' @export
#' @importFrom tinyplot tinyplot type_text
plot_mov <- function(fit, s = 1, y, a, palette = "Peach") {

  dat <- get_MSAdata(fit)

  nm <- dat@Dmodel@nm
  nr <- dat@Dmodel@nr
  if (missing(y)) y <- dat@Dmodel@ny
  if (missing(a)) a <- 1
  rname <- dat@Dlabel@region
  mname <- dat@Dlabel@season

  mov <- array(fit@report$mov_ymarrs[y, , a+1, , , s], c(nm, nr, nr))

  #if (nm > 1) {
  #  old_mar <- par()$mar
  #  old_mfrow = par()$mfrow
  #  par(mar = c(4, 4, 1, 1))
  #  on.exit(par(mar = old_mar, mfrow = old_mfrow))
  #  par(mfrow = c(2, ceiling(nm/2)))
  #}

  dist_eq <- calc_eqdist(mov, start = fit@report$recdist_rs[, s], m_start = dat@Dstock@m_spawn)

  df_mov <- structure(mov, dimnames = list(Season = mname, Origin = 1:nr, Destination = 1:nr)) %>%
    reshape2::melt()

  df_eq <- structure(dist_eq, dimnames = list(Season = mname, Origin = 1:nr)) %>%
    reshape2::melt() %>%
    cbind("Destination" = nr + 1.5)

  df <- rbind(df_mov, df_eq[, c("Season", "Origin", "Destination", "value")])
  #df$label <- format(round(df$value, 2), nsmall = 2)
  df$label <- round(df$value, 2)

  tick_fn <- function(i) ifelse(i > nr, "Eq.", as.character(rname[i]))
  tinyplot_args <- list(
    xmin = df$Destination - 0.5, xmax = df$Destination + 0.5, ymin = df$Origin - 0.5, ymax = df$Origin + 0.5,
    by = df$value,
    facet = df$Season, xlab = "Destination", ylab = "Origin",
    bg = "by", col = "black",
    legend = substitute(legend(title = "Proportion")),
    yaxs = "i", xaxs = "i",
    yaxl = tick_fn, xaxl = tick_fn,
    yaxb = 1:nr, xaxb = c(1:nr, nr + 1.5),
    type = "rect", palette = palette
  )
  do.call(tinyplot, tinyplot_args)

  tinyplot(
    x = df$Destination, y = df$Origin, facet = df$Season,
    type = type_text(labels = df$label, adj = 0.5),
    add = TRUE
  )

  #for(m in 1:nm) {
  #  .plot_mov(m = mov[m, , ], p = dist_eq[m, ], rname = rname, xlab = "", ylab = "", palette = palette)
  #  if (nm > 1) title(mname[m], font.main = 1)
  #}
  #par(mfrow = c(1, 1))
  #mtext("Destination", side = 1, line = 3.5)
  #mtext("Origin", side = 2, line = 3)

  invisible()
}

#' @rdname plot-MSA-state
#' @aliases plot_recdist
#' @details
#' - `plot_recdist` plots the distribution of recruitment for each stock
#' @export
plot_recdist <- function(fit, palette = "Peach") {
  dat <- get_MSAdata(fit)

  nr <- dat@Dmodel@nr

  if (nr > 1) {
    ns <- dat@Dmodel@ns
    rname <- dat@Dlabel@region
    sname <- dat@Dlabel@stock

    recdist <- fit@report$recdist_rs

    #vcol <- hcl.colors(100, palette)
#
    #graphics::plot.default(
    #  NULL, xlab = "Stock", ylab = "Region", xaxs = "i", yaxs = "i",
    #  xaxt = "n", yaxt = "n", xlim = c(1, ns+1), ylim = c(1, nr+1)
    #)
    #for(x in 1:ns) {
    #  for(y in 1:nr) {
    #    m_yx <- round(recdist[y, x], 2)
    #    rect(xleft = x, ybottom = y, xright = x+1, ytop = y+1, col = vcol[100 * m_yx])
    #    text(x + 0.5, y + 0.5, m_yx)
    #  }
    #}
#
    #axis(1, at = 1:ns + 0.5, labels = as.character(sname), font = 2, cex.axis = 0.75)
    #axis(2, at = 1:nr + 0.5, labels = as.character(rname), font = 2, cex.axis = 0.75)

    df <- structure(recdist, dimnames = list(r = 1:nr, s = 1:ns)) %>%
      reshape2::melt()
    df$label <- round(df$value, 2)

    tick_fnx <- function(i) sname[i]
    tick_fny <- function(i) rname[i]

    tinyplot_args <- list(
      xmin = df$s - 0.5, xmax = df$s + 0.5, ymin = df$r - 0.5, ymax = df$r + 0.5,
      by = df$value, xlab = "Stock", ylab = "Region",
      bg = "by", col = "black", border = "black",
      legend = substitute(legend(title = "Proportion")),
      yaxs = "i", xaxs = "i",
      yaxl = tick_fny, xaxl = tick_fnx,
      yaxb = 1:nr, xaxb = 1:ns,
      type = "rect", palette = palette
    )
    do.call(tinyplot, tinyplot_args)

    tinyplot(
      x = df$s, y = df$r,
      type = type_text(labels = df$label, adj = 0.5),
      add = TRUE
    )
  }

  invisible()
}

# #' @importFrom grDevices hcl.colors
# #' @importFrom graphics rect text
# .plot_mov <- function(m, p, xlab = "Destination", ylab = "Origin",
#                       nr = length(p), rname = paste("Region", 1:nr), palette = "Peach") {
#
#   vcol <- hcl.colors(100, palette)
#
#   graphics::plot.default(
#     NULL, xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "i",
#     xaxt = "n", yaxt = "n", xlim = c(1, nr+3), ylim = c(1, nr+1)
#   )
#   for(x in 1:nr) {
#     for(y in 1:nr) {
#       m_yx <- round(m[y, x], 2)
#       rect(xleft = x, ybottom = y, xright = x+1, ytop = y+1, col = vcol[100 * m_yx])
#       text(x + 0.5, y + 0.5, m_yx)
#     }
#   }
#   eq_val <- round(p, 2)
#   eq_col <- rep(NA, length(p))
#   eq_col[eq_val > 0] <- vcol[100 * eq_val]
#   rect(nr + 2, ybottom = 1:nr, xright = nr + 3, ytop = 1:nr + 1, col = eq_col)
#   text(nr + 2.5, 1:nr + 0.5, eq_val)
#
#   axis(2, at = 1:nr + 0.5, labels = as.character(rname), font = 2, cex.axis = 0.75)
#   axis(1, at = c(1:nr, nr+2) + 0.5, labels = c(as.character(rname), "Eq."), font = 2, cex.axis = 0.75)
#
#   invisible()
# }
