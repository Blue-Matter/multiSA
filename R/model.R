

#' Fit MSA model
#'
#' Wrapper function that calls RTMB to create the model and perform the numerical optimization
#'
#' @param MSAdata Data object. Class [MSAdata-class], validated by [check_data()]
#' @param parameters List of parameters, e.g., returned by [make_parameters()] and validated by [check_parameters()].
#' @param map List of parameters indicated whether they are fixed and how they are shared, e.g., returned by [make_parameters()].
#' See [RTMB::MakeADFun()].
#' @param random Character vector indicating the parameters that are random effects, e.g., returned by [make_parameters()].
#' @param run_model Logical, whether to fit the model through [stats::nlminb()].
#' @param do_sd Logical, whether to calculate the standard errors with [RTMB::sdreport()].
#' @param report Logical, whether to return the report list with `obj$report(obj$env$last.par.best)`.
#' @param silent Logical, whether to report progress to console. **Not passed to [TMB::MakeADFun()].**
#' @param control Passed to [stats::nlminb()]
#' @param ... Other arguments to [RTMB::MakeADFun()].
#' @returns A [MSAassess-class] object.
#' @importFrom methods new
#' @seealso [report()] [retrospective()]
#' @export
fit_MSA <- function(MSAdata, parameters, map = list(), random = NULL,
                     run_model = TRUE, do_sd = TRUE, report = TRUE, silent = FALSE,
                     control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {

  MSAdata@Misc$map <- map
  MSAdata@Misc$random <- random

  old_comparison <- TapeConfig()["comparison"]

  on.exit(TapeConfig(comparison = old_comparison))
  TapeConfig(comparison = "tape")

  func <- function(p) .MSA(p, d = MSAdata)

  if (!silent) message("Building model..")
  obj <- RTMB::MakeADFun(
    func = func, parameters = parameters,
    map = map, random = random,
    silent = TRUE,
    ...
  )

  if (!silent) {
    fn <- obj$fn()
    if (is.na(fn)) {
      message_oops("Objective function is NA at initial values.")
      report_start <- obj$report()
      if (MSAdata@Dmodel@condition == "catch" && any(is.na(report_start$F_ymfr))) {
        message_oops("NA's found in F array. Try increasing start value of R0.")
      }

    } else if (is.infinite(fn)) {
      message_oops("Objective function is infinite at initial values.")

    }

    gr <- obj$gr()
    if (any(is.na(gr))) message_oops("Gradients of NA at initial values were found.")
    if (any(!gr, na.rm = TRUE)) {
      message_oops("Gradients of zero at initial values, can be indicative of over-parameterization or non-identifiable parameters.")
    }
  }

  M <- new("MSAassess", obj = obj)

  if (run_model) {
    m <- optimize_RTMB(obj, do_sd = do_sd, control = control, silent = silent)
    if (is.character(m$opt)) {
      message_oops("Error message from optimization:\n", m$opt)
    } else {
      M@opt <- m$opt
    }
    if (do_sd) M@SD <- m$SD
  }

  if (report) {
    if (!silent) message("Generating report list..")
    M@report <- obj$report(obj$env$last.par.best) %>% update_report(MSAdata)
  }
  if (!silent) message("Complete.")
  return(M)
}

update_report <- function(r, MSAdata) {

  if (is.null(r$F_yas)) {
    nr <- MSAdata@Dmodel@nr
    nm <- MSAdata@Dmodel@nm

    ny <- MSAdata@Dmodel@ny
    na <- MSAdata@Dmodel@na
    ns <- MSAdata@Dmodel@ns

    if (nr == 1 && nm == 1) {
      r$F_yas <- array(r$F_ymars[, 1, , 1, ], c(ny, na, ns))
    } else {
      Fmax <- MSAdata@Dmodel@Fmax
      nf <- MSAdata@Dfishery@nf

      r$F_yas <- sapply2(1:ns, function(s) {
        sapply(1:na, function(a) {
          sapply(1:ny, function(y) {
            calc_summary_F(M = r$M_yas[y, a, s], N = sum(r$N_ymars[y, 1, a, , s]),
                           CN = sum(r$CN_ymafrs[y, , a, , , s]), Fmax = nf * Fmax)
          })
        })
      })

      if (any(is.na(r$F_yas))) message_oops("NA in predicted catch at age")
      if (any(is.infinite(r$F_yas))) message_oops("Catch-at-age exceeds abundance-at-age")
    }
    r$Z_yas <- r$F_yas + r$M_yas
  }
  return(r)
}

.MSA <- function(p = list(), d) {

  # Dispatch method for AD variables ----
  is_ad <- any(sapply(p, inherits, "advector"))
  if (is_ad) {
    `[<-` <- RTMB::ADoverload("[<-")
  }

  # Assign data variables to environment, see OBS() for simulation ----
  #getAllS4(d@Dmodel, d@Dstock, d@Dfishery, d@Dsurvey, d@DCKMR, d@Dtag)
  Dmodel <- d@Dmodel
  Dstock <- d@Dstock
  Dfishery <- d@Dfishery
  Dsurvey <- d@Dsurvey
  DCKMR <- d@DCKMR
  Dtag <- d@Dtag

  ny <- Dmodel@ny
  nm <- Dmodel@nm
  na <- Dmodel@na
  nl <- Dmodel@nl
  nr <- Dmodel@nr
  ns <- Dmodel@ns

  nf <- Dfishery@nf
  ni <- Dsurvey@ni

  map <- d@Misc$map
  random <- d@Misc$random

  # Population arrays ----
  N_ymars <- array(NA_real_, c(ny + 1, nm, na, nr, ns))
  B_ymrs <- array(NA_real_, c(ny, nm, nr, ns))

  Npsp_yars <-
    Nsp_yars <- array(NA_real_, c(ny, na, nr, ns))
  S_yrs <- array(NA_real_, c(ny, nr, ns))
  R_ys <-
    Rdev_ys <- array(NA_real_, c(ny, ns))
  F_ymars <-
    Z_ymars <- array(NA_real_, c(ny, nm, na, nr, ns))
  M_yas <- mat_yas <- array(NA_real_, c(ny, na, ns))
  if (length(DCKMR@HSP_s)) F_yas <- array(NA_real_, c(ny, na, ns))
  mov_ymarrs <- array(NA_real_, c(ny, nm, na, nr, nr, ns))
  recdist_rs <- array(NA_real_, c(nr, ns))

  # Fishery arrays ----
  sel_ymafs <- array(NA_real_, c(ny, nm, na, nf, ns))

  F_ymfr <-
    log_F_ymfr <- array(NA_real_, c(ny, nm, nf, nr))
  F_ymafrs <-
    CN_ymafrs <- array(NA_real_, c(ny, nm, na, nf, nr, ns))
  if (any(Dfishery@CALobs_ymlfr > 0, na.rm = TRUE)) CN_ymlfrs <- array(NA_real_, c(ny, nm, nl, nf, nr, ns))
  CB_ymfrs <-
    VB_ymfrs <- array(NA_real_, c(ny, nm, nf, nr, ns))

  # Index of abundance arrays ----
  if (ni > 0) {
    IN_ymais <-
      sel_ymais <- array(NA_real_, c(ny, nm, na, ni, ns))
    if (any(Dsurvey@IALobs_ymli > 0, na.rm = TRUE)) IN_ymlis <- array(NA_real_, c(ny, nm, nl, ni, ns))
    VI_ymi <- I_ymi <- array(NA_real_, c(ny, nm, ni))
  }

  # Transform parameters ----
  ## Maturity at age ogive ----
  if (is.null(map$mat_ps)) map$mat_ps <- matrix(TRUE, 2, ns)
  map$mat_ps <- matrix(as.character(map$mat_ps), 2, ns)
  mat_yas[] <- sapply2(1:ns, function(s) {
    if (all(is.na(map$mat_ps[, s]))) {
      Dstock@mat_yas[1:ny, , s]
    } else {
      m <- conv_mat(p$mat_ps[, s], na)
      matrix(m, ny, na, byrow = TRUE)
    }
  })

  ## Natural mortality ----
  for(s in 1:ns) {
    if (!is.null(map$log_M_s) && is.na(map$log_M_s[s])) {
      M_yas[, , s] <- Dstock@M_yas[1:ny, , s]
    } else {
      M_yas[, , s] <- matrix(exp(p$log_M_s[s]), ny, na)
    }
  }

  ## Fishery selectivity ----
  q_fs <- exp(p$log_q_fs)
  selconv_pf <- conv_selpar(p$sel_pf, type = Dfishery@sel_f, maxage = na - 1, maxL = 0.95 * max(Dmodel@lmid))
  sel_lf <- calc_sel_len(selconv_pf, Dmodel@lmid, type = Dfishery@sel_f)

  ## Fishing mortality ----
  if (Dmodel@condition == "F") {
    log_Fmult_f <- sapply(1:nf, function(f) p$log_Fdev_ymfr[Dmodel@y_Fmult_f[f], Dmodel@m_Fmult_f[f], f, Dmodel@r_Fmult_f[f]])
    log_F_ymfr[] <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply(1:nm, function(m) {
          sapply(1:ny, function(y) {
            Fmult_y <- y == Dmodel@y_Fmult_f[f]
            Fmult_m <- m == Dmodel@m_Fmult_f[f]
            Fmult_r <- r == Dmodel@r_Fmult_f[f]
            if (Fmult_y && Fmult_m && Fmult_r) {
              log_Fmult_f[f]
            } else {
              log_Fmult_f[f] + p$log_Fdev_ymfr[y, m, f, r]
            }
          })
        })
      })
    })
    F_ymfr[] <- exp(log_F_ymfr)
  }

  ## Index selectivity ----
  if (ni > 0) {
    selconv_pi <- conv_selpar(p$sel_pi, type = Dsurvey@sel_i, maxage = na - 1, maxL = 0.95 * max(Dmodel@lmid))
    sel_li <- calc_sel_len(selconv_pi, Dmodel@lmid, type = Dsurvey@sel_i)
  }

  ## Fishery and index selectivity ----
  # Check for fishery selectivity blocks
  tv_flensel <- any(
    sapply(1:nf, function(f) {
      bb <- Dfishery@sel_block_yf[, f]
      lensel <- any(grepl("length", Dfishery@sel_f[bb]))
      change_sel <- length(unique(bb)) > 1
      change_sel || lensel
    })
  )
  if (ni > 0) { # Check for length selectivity blocks
    ilensel <- any(sapply(1:ni, function(i) grepl("length", Dsurvey@sel_i[i])))
  }

  for (s in 1:ns) {
    # Fishery selectivity

    # Check for time-varying growth
    if (tv_flensel) {
      tv_growth <- any(sapply(2:ny, function(y) max(Dstock@LAK_ymals[y, , , , s] - Dstock@LAK_ymals[1, , , , s])) > 0)
      tv_fagesel_growth <- tv_growth
    } else {
      tv_fagesel_growth <- FALSE
    }

    # Check for time-varying maturity
    tv_mat <- any(sapply(2:ny, function(y) max(mat_yas[y, , s] - mat_yas[1, , s])) > 0)
    tv_fsel_mat <- tv_mat && any(Dfishery@sel_f == "SB")

    if (tv_fagesel_growth || tv_fsel_mat) {
      for (y in 1:ny) {
        for (m in 1:nm) {
          sel_ymafs[y, m, , , s] <- calc_fsel_age(
            sel_lf, Dstock@LAK_ymals[y, m, , , s], Dfishery@sel_f, selconv_pf, Dfishery@sel_block_yf[y, ], mat = mat_yas[y, , s], a = seq(1, na) - 1
          )
        }
      }
    } else {
      for (m in 1:nm) {
        sel_ymafs[1, m, , , s] <- calc_fsel_age(
          sel_lf, Dstock@LAK_ymals[1, m, , , s], Dfishery@sel_f, selconv_pf, Dfishery@sel_block_yf[1, ], mat = mat_yas[1, , s], a = seq(1, na) - 1
        )
      }
      fsel_ind <- fsel1_ind <- as.matrix(expand.grid(y = 2:ny, m = 1:m, a = 1:na, f = 1:nf, s = s))
      fsel1_ind[, "y"] <- 1
      sel_ymafs[fsel_ind] <- sel_ymafs[fsel1_ind]
    }

    # Survey selectivity
    if (ni > 0) {
      tv_iagesel_growth <- ilensel && tv_growth
      tv_isel_mat <- tv_mat && any(Dsurvey@sel_i == "SB")
      if (tv_iagesel_growth || tv_isel_mat) {
        for (y in 1:ny) {
          for (m in 1:nm) {
            sel_ymais[y, m, , , s] <- calc_isel_age(
              sel_li, Dstock@LAK_ymals[y, m, , , s], Dsurvey@sel_i, selconv_pi, sel_ymafs[y, m, , , s], mat = mat_yas[y, , s], a = seq(1, na) - 1
            )
          }
        }
      } else {
        for (m in 1:nm) {
          sel_ymais[1, m, , , s] <- calc_isel_age(
            sel_li, Dstock@LAK_ymals[1, m, , , s], Dsurvey@sel_i, selconv_pi, sel_ymafs[1, m, , , s], mat = mat_yas[1, , s], a = seq(1, na) - 1
          )
        }
        isel_ind <- isel1_ind <- as.matrix(expand.grid(y = 2:ny, m = 1:m, a = 1:na, i = 1:ni, s = s))
        isel1_ind[, "y"] <- 1
        sel_ymais[isel_ind] <- sel_ymais[isel1_ind]
      }
    }
  }

  ## Stock distribution and movement parameters ----
  recdist_rs[] <- sapply(1:ns, function(s) softmax(p$log_recdist_rs[, s]))

  # Movement
  if (nr > 1) {
    for (yy in 1:nrow(Dtag@tag_yy)) {
      yvec <- which(Dtag@tag_yy[yy, ] > 0)
      y1 <- yvec[1]

      for (m in 1:nm) {
        mov_ymarrs[y1, m, , , , ] <- sapply2(1:ns, function(s) {
          mov_arr <- array(0, c(na, nr, nr))
          r_eff <- Dstock@presence_rs[, s]
          nr_eff <- sum(Dstock@presence_rs[, s])
          mov_arr[, r_eff, r_eff] <- conv_mov(
            p$mov_x_marrs[m, , r_eff, r_eff, s], p$mov_g_ymars[y1, m, , r_eff, s], p$mov_v_ymas[y1, m, , s], na, nr_eff
          )
          return(mov_arr)
        })
      }

      if (length(yvec) > 1) {
        mov_ind <- mov1_ind <- as.matrix(expand.grid(y = yvec[-1], m = 1:m, a = 1:na, rf = 1:nr, rt = 1:nr, s = 1:ns))
        mov1_ind[, "y"] <- yvec[1]
        mov_ymarrs[mov_ind] <- mov_ymarrs[mov1_ind]
      }
    }
  } else {
    mov_ymarrs[] <- 1
  }

  ## Stock recruit parameters ----
  R0_s <- exp(p$t_R0_s) * Dmodel@scale_s
  h_s <- sapply(1:ns, function(s) conv_steepness(p$t_h_s[s], Dstock@SRR_s[s]))
  kappa_s <- sapply(1:ns, function(s) SRkconv(h_s[s], Dstock@SRR_s[s]))

  if (nr == 1 && nm == 1) {
    nyinit <- 1L
    initNPR0_yars <- sapply2(1:ns, function(s) calc_NPR(M_yas[Dmodel@y_phi, , s])) %>% array(c(nyinit, na, nr, ns))
    phi_s <- sapply(1:ns, function(s) {
      calc_phi_simple(M_yas[Dmodel@y_phi, , s], mat_a = mat_yas[Dmodel@y_phi, , s], fec_a = Dstock@fec_yas[Dmodel@y_phi, , s],
                      delta = Dstock@delta_s[s])
    })
  } else {
    nyinit <- Dmodel@nyinit
    NPR_unfished <- calc_phi_project(
      nyinit, nm, na, nf = 1, nr, ns, M_as = M_yas[Dmodel@y_phi, , ], mov_marrs = mov_ymarrs[Dmodel@y_phi, , , , , ],
      mat_as = mat_yas[Dmodel@y_phi, , ], fec_as = Dstock@fec_yas[Dmodel@y_phi, , ], m_spawn = Dstock@m_spawn, m_advanceage = Dstock@m_advanceage,
      delta_s = Dstock@delta_s, natal_rs = Dstock@natal_rs, recdist_rs = recdist_rs
    )
    initNPR0_yars <- array(NPR_unfished[["N_ymars"]][1:nyinit, 1, , , ], c(nyinit, na, nr, ns))
    phi_s <- sapply(1:ns, function(s) sum(NPR_unfished[["S_yrs"]][nyinit, , s]))
  }

  SB0_s <- R0_s * phi_s
  sralpha_s <- kappa_s/phi_s
  srbeta_s <- sapply(1:ns, function(s) SRbetaconv(h_s[s], R0_s[s], phi_s[s], SRR = Dstock@SRR_s[s]))

  ## Recruitment deviates ----
  sdr_s <- exp(p$log_sdr_s)
  bcr_s <- -0.5 * sdr_s * sdr_s

  par_rdev_ys <- matrix(TRUE, ny, ns)
  bcr_ys <- sapply(1:ns, function(s) bcr_s[s] * Dmodel@pbc_rdev_ys[, s])
  if (!is.null(map$log_rdev_ys)) {
    par_rdev_ys[] <- matrix(!is.na(map$log_rdev_ys) & !duplicated(map$log_rdev_ys, MARGIN = 0), ny, ns)
    bcr_ys[is.na(map$log_rdev_ys)] <- 0
  }

  par_initrdev_as <- matrix(TRUE, na-1, ns)
  bcrinit_as <- sapply(1:ns, function(s) bcr_s[s] * Dmodel@pbc_initrdev_as[, s])
  if (!is.null(map$log_initrdev_as)) {
    par_initrdev_as[] <- matrix(!is.na(map$log_initrdev_as) & !duplicated(map$log_initrdev_as, MARGIN = 0), na-1, ns)
    bcrinit_as[is.na(map$log_initrdev_as)] <- 0
  }
  Rdev_ys[] <- exp(p$log_rdev_ys + bcr_ys)

  ## Miscellaneous penalty term, e.g., F > Fmax
  penalty <- 0

  # First year, first season initialization ----
  initRdev_as <- matrix(NA_real_, na, ns)
  initRdev_as[-1, ] <- exp(p$log_initrdev_as + bcrinit_as)
  initRdev_as[1, ] <- Rdev_ys[1, ]
  initF_mfr <- exp(p$log_initF_mfr)

  initNPR_yars <- array(NA_real_, c(nyinit, na, nr, ns))
  initN_ars <- array(NA_real_, c(na, nr, ns))

  initZ_mars <- array(NA_real_, c(nm, na, nr, ns))

  initCN_mafrs <- array(NA_real_, c(nm, na, nf, nr, ns))
  initCB_mfrs <- array(NA_real_, c(nm, nf, nr, ns))

  if (all(Dfishery@Cinit_mfr < 1e-8)) {
    initNPR_yars[] <- initNPR0_yars
    initphi_s <- phi_s
    initR_s <- R0_s
  } else {
    NPR_init <- calc_phi_project( #nyinit = 1 if nm == 1 && nr == 1
      nyinit, nm, na, nf, nr, ns, F_mfr = initF_mfr, sel_mafs = sel_ymafs[1, , , , ],
      fwt_mafs = Dfishery@fwt_ymafs[1, , , , ], q_fs = q_fs,
      M_as = M_yas[1, , ], mov_marrs = mov_ymarrs[Dmodel@y_phi, , , , , ],
      mat_as = mat_yas[1, , ], fec_as = Dstock@fec_yas[1, , ], m_spawn = Dstock@m_spawn, m_advanceage = Dstock@m_advanceage,
      delta_s = Dstock@delta_s, natal_rs = Dstock@natal_rs, recdist_rs = recdist_rs
    )
    initNPR_yars[] <- NPR_init[["N_ymars"]][1:nyinit, 1, , , ]
    initphi_s <- sapply(1:ns, function(s) sum(NPR_init[["S_yrs"]][nyinit, , s]))
    initZ_mars[] <- NPR_init[["Z_ymars"]][nyinit, , , , ]

    initCN_mafrs[] <- NPR_init[["CN_ymafrs"]][nyinit, , , , , ]
    initCB_mfrs[] <- NPR_init[["CB_ymfrs"]][nyinit, , , , ]

    initR_s <- sapply(1:ns, function(s) {
      calc_recruitment(initphi_s[s], Dstock@SRR_s[s], eq = TRUE, a = sralpha_s[s], b = srbeta_s[s])
    })
  }

  initN_ars[] <- sapply2(1:ns, function(s) initR_s[s] * initRdev_as[, s] * initNPR_yars[nyinit, , , s])

  # Run population model ----
  pop <- calc_population(
    ny, nm, na, nf, nr, ns, initN_ars, mov_ymarrs, M_yas, Dstock@SRR_s, sralpha_s, srbeta_s,
    mat_yas, Dstock@fec_yas, Rdev_ys, Dstock@m_advanceage, Dstock@m_spawn, Dstock@delta_s, Dstock@natal_rs, recdist_rs = recdist_rs,
    Dfishery@fwt_ymafs, q_fs, sel_ymafs,
    condition = Dmodel@condition, F_ymfr = F_ymfr, Cobs_ymfr = Dfishery@Cobs_ymfr, Fmax = Dmodel@Fmax, nitF = Dmodel@nitF
  )

  # Assign population arrays ----
  N_ymars[] <- pop$N_ymars
  F_ymars[] <- pop$F_ymars
  Z_ymars[] <- pop$Z_ymars
  F_ymafrs[] <- pop$F_ymafrs
  CN_ymafrs[] <- pop$CN_ymafrs
  CB_ymfrs[] <- pop$CB_ymfrs
  VB_ymfrs[] <- pop$VB_ymfrs
  Nsp_yars[] <- pop$Nsp_yars
  Npsp_yars[] <- pop$Npsp_yars
  S_yrs[] <- pop$S_yrs
  R_ys[] <- pop$R_ys
  penalty <- penalty + pop$penalty

  if (Dmodel@condition == "catch") F_ymfr[] <- pop$F_ymfr

  ind_ymars <- as.matrix(expand.grid(y = 1:ny, m = 1:nm, a = 1:na, r = 1:nr, s = 1:ns))
  ymas_ymars <- ind_ymars[, c("y", "m", "a", "s")]

  B_ymrs[] <- array(N_ymars[ind_ymars] * Dstock@swt_ymas[ymas_ymars], c(ny, nm, na, nr, ns)) %>%
    apply(c(1, 2, 4, 5), sum)

  # Likelihoods ----
  y_like <- seq(1, ny - Dmodel@nyret)

  ## Initial catch ----
  Cinit_mfr <- Dfishery@Cinit_mfr
  any_Cinit <- any(Cinit_mfr >= 1e-8)
  if (any_Cinit) {
    initCB_mfr <- apply(initCB_mfrs, 1:3, sum)

    Cinit_mfr <- OBS(Cinit_mfr)
    loglike_Cinit_mfr <- dnorm(log(Cinit_mfr/initCB_mfr), 0, 0.01, log = TRUE)
    loglike_Cinit_mfr[Cinit_mfr < 1e-8] < 0
  } else {
    loglike_Cinit_mfr <- 0
  }

  ## Catch ----
  Cobs_ymfr <- Dfishery@Cobs_ymfr
  if (Dmodel@condition == "F") {
    CB_ymfr <- apply(CB_ymfrs, 1:4, sum)
    Cobs_ymfr <- OBS(Cobs_ymfr)

    loglike_Cobs_ymfr <- array(0, c(ny, nm, nf, nr))
    loglike_Cobs_ymfr[] <- dnorm(log(Cobs_ymfr/CB_ymfr), 0, Dfishery@Csd_ymfr, log = TRUE)
    loglike_Cobs_ymfr[Cobs_ymfr < 1e-8] <- 0
    loglike_Cobs_ymfr[1:ny > max(y_like), , , ] <- 0
  } else {
    loglike_Cobs_ymfr <- 0
  }

  ## Marginal fishery age composition ----
  CAAobs_ymafr <- Dfishery@CAAobs_ymafr
  any_CAA <- any(CAAobs_ymafr > 0, na.rm = TRUE)
  if (any_CAA) {
    CN_ymafr <- apply(CN_ymafrs, 1:5, sum)
    CAAobs_ymafr <- OBS(CAAobs_ymafr)

    loglike_CAA_ymfr <- array(0, c(ny, nm, nf, nr))
    loglike_CAA_ymfr[y_like, , , ] <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply(1:nm, function(m) {
          sapply(y_like, function(y) {
            pred <- CN_ymafr[y, m, , f, r]
            like_comp(obs = (Cobs_ymfr[y, m, f, r] >= 1e-8) * CAAobs_ymafr[y, m, , f, r],
                      pred = pred, type = Dfishery@fcomp_like,
                      N = Dfishery@CAAN_ymfr[y, m, f, r], theta = Dfishery@CAAtheta_f[f],
                      stdev = sqrt(sum(pred)/pred))
          })
        })
      })
    })
  } else {
    loglike_CAA_ymfr <- 0
  }

  ## Marginal fishery length composition ----
  CALobs_ymlfr <- Dfishery@CALobs_ymlfr
  any_CAL <- any(CALobs_ymlfr > 0, na.rm = TRUE)
  if (any_CAL) {
    for(y in 1:ny) {
      for(m in 1:nm) {
        CN_ymlfrs[y, m, , , , ] <- sapply2(1:ns, function(s) {
          sapply2(1:nr, function(r) {
            sapply(1:nf, function(f) CN_ymafrs[y, m, , f, r, s] %*% Dstock@LAK_ymals[y, m, , , s])
          })
        })
      }
    }
    CN_ymlfr <- apply(CN_ymlfrs, 1:5, sum)
    CALobs_ymlfr <- OBS(CALobs_ymlfr)

    loglike_CAL_ymfr <- array(0, c(ny, nm, nf, nr))
    loglike_CAL_ymfr[y_like, , , ] <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply(1:nm, function(m) {
          sapply(y_like, function(y) {
            pred <- CN_ymlfr[y, m, , f, r]
            like_comp(obs = (Cobs_ymfr[y, m, f, r] >= 1e-8) * CALobs_ymlfr[y, m, , f, r],
                      pred = pred, type = Dfishery@fcomp_like,
                      N = Dfishery@CALN_ymfr[y, m, f, r], theta = Dfishery@CALtheta_f[f],
                      stdev = sqrt(sum(pred)/pred))
          })
        })
      })
    })
  } else {
    loglike_CAL_ymfr <- 0
  }

  ## Index ----
  Iobs_ymi <- Dsurvey@Iobs_ymi
  if (ni > 0) {
    for(y in 1:ny) {
      for(m in 1:nm) {
        IN_ymais[y, m, , , ] <- calc_index(
          N = N_ymars[y, m, , , ], Z = Z_ymars[y, m, , , ], sel = sel_ymais[y, m, , , ],
          na = na, nr = nr, ns = ns, ni = ni, samp = Dsurvey@samp_irs, delta = Dsurvey@delta_i
        )
        VI_ymi[y, m, ] <- sapply(1:ni, function(i) {
          I_s <- sapply(1:ns, function(s) {
            w <- if (Dsurvey@unit_i[i] == "N") 1 else Dstock@swt_ymas[y, m, , s]
            sum(IN_ymais[y, m, , i, s] * w)
          })
          sum(I_s)
        })
      }
    }
    q_i <- sapply(1:ni, function(i) calc_q(Iobs_ymi[, , i], B = VI_ymi[, , i]))
    I_ymi[] <- sapply2(1:ni, function(i) q_i[i] * VI_ymi[, , i])

    Iobs_ymi <- OBS(Iobs_ymi)
    loglike_I_ymi <- array(0, c(ny, nm, ni))
    loglike_I_ymi[] <- dnorm(log(Iobs_ymi/I_ymi), 0, Dsurvey@Isd_ymi, log = TRUE)
    loglike_I_ymi[is.na(Iobs_ymi)] <- 0
    loglike_I_ymi[1:ny > max(y_like), , ] <- 0

  } else {
    loglike_I_ymi <- 0
  }

  IAAobs_ymai <- Dsurvey@IAAobs_ymai
  any_IAA <- ni > 0 && any(IAAobs_ymai > 0, na.rm = TRUE)
  if (any_IAA) {
    IN_ymai <- apply(IN_ymais, 1:4, sum)
    IAAobs_ymai <- OBS(IAAobs_ymai)

    loglike_IAA_ymi <- array(0, c(ny, nm, ni))
    loglike_IAA_ymi[y_like, , ] <- sapply2(1:ni, function(i) {
      sapply(1:nm, function(m) {
        sapply(y_like, function(y) {
          pred <- IN_ymai[y, m, , i]
          like_comp(obs = IAAobs_ymai[y, m, , i], pred = pred, type = Dsurvey@icomp_like,
                    N = Dsurvey@IAAN_ymi[y, m, i], theta = Dsurvey@IAAtheta_i[i],
                    stdev = sqrt(sum(pred)/pred))
        })
      })
    })
  } else {
    loglike_IAA_ymi <- 0
  }

  IALobs_ymli <- Dsurvey@IALobs_ymli
  any_IAL <- ni > 0 && any(IALobs_ymli > 0, na.rm = TRUE)
  if (any_IAL) {
    for(y in 1:ny) {
      for(m in 1:nm) {
        IN_ymlis[y, m, , , ] <- sapply2(1:ns, function(s) {
          sapply(1:ni, function(i) IN_ymais[y, m, , i, s] %*% Dstock@LAK_ymals[y, m, , , s])
        })
      }
    }
    IN_ymli <- apply(IN_ymlis, 1:4, sum)
    IALobs_ymli <- OBS(IALobs_ymli)

    loglike_IAL_ymi <- array(0, c(ny, nm, ni))
    loglike_IAL_ymi[y_like, , ] <- sapply2(1:ni, function(i) {
      sapply(1:nm, function(m) {
        sapply(y_like, function(y) {
          pred <- IN_ymli[y, m, , i]
          like_comp(obs = IALobs_ymli[y, m, , i], pred = pred, type = Dsurvey@icomp_like,
                    N = Dsurvey@IALN_ymi[y, m, i], theta = Dsurvey@IALtheta_i[i],
                    stdev = sqrt(sum(pred)/pred))
        })
      })
    })
  } else {
    loglike_IAL_ymi <- 0
  }

  ## Stock composition ----
  SC_ymafrs <- Dfishery@SC_ymafrs
  any_SC <- ns > 1 && length(SC_ymafrs)
  if (any_SC) {
    SC_ymafrs <- OBS(SC_ymafrs)
    SCpred_ymafrs <- sapply2(1:nrow(Dfishery@SC_ff), function(ff) {
      fvec <- Dfishery@SC_ff[ff, ]
      sapply2(1:nrow(Dfishery@SC_aa), function(aa) {
        avec <- Dfishery@SC_aa[aa, ]
        apply(CN_ymafrs[, , avec, fvec, , , drop = FALSE], c(1, 2, 5, 6), sum)
      })
    }) %>%
      aperm(c(1, 2, 5, 6, 3, 4))

    loglike_SC_ymafr <- array(0, dim(SC_ymafrs)[1:5])
    loglike_SC_ymafr[y_like, , , , ] <- sapply2(1:nr, function(r) {
      sapply2(1:nrow(Dfishery@SC_ff), function(ff) {
        sapply2(1:nrow(Dfishery@SC_aa), function(aa) {
          sapply(1:nm, function(m) {
            sapply(y_like, function(y) {
              pred <- SCpred_ymafrs[y, m, aa, ff, r, ]
              Cobs <- sum(Cobs_ymfr[y, m, ff, r])
              like_comp(obs = (Cobs >= 1e-8) * SC_ymafrs[y, m, aa, ff, r, ],
                        pred = CondExpGt(pred, 1e-8, pred, 1e-8), type = Dfishery@SC_like,
                        N = Dfishery@SCN_ymafr[y, m, aa, ff, r], theta = Dfishery@SCtheta_f[ff],
                        stdev = Dfishery@SCstdev_ymafrs[y, m, aa, ff, r, ])
            })
          })
        })
      })
    })
  } else {
    loglike_SC_ymafr <- 0
  }

  ## CKMR ----
  if (length(DCKMR@POP_s)) {
    pPOP_s <- lapply(1:ns, function(s) {
      calc_POP(t = DCKMR@POP_s[[s]]$t, a = DCKMR@POP_s[[s]]$a, y = DCKMR@POP_s[[s]]$y,
               N = apply(Nsp_yars[, , , s], 1:2, sum), fec = Dstock@fec_yas[, , s])
    })
    loglike_POP_s <- lapply(1:ns, function(s) {
      val <- 0
      if (DCKMR@POP_s[[s]]$y %in% y_like) {
        val <- like_CKMR(n = DCKMR@POP_s[[s]]$n, m = DCKMR@POP_s[[s]]$m, p = pPOP_s[[s]], type = DCKMR@CKMR_like)
      }
      return(val)
    })
  } else {
    loglike_POP_s <- 0
  }

  if (length(DCKMR@HSP_s)) {
    ## Summary F and Z by year ----
    F_yas[] <- sapply2(1:ns, function(s) {
      sapply(1:na, function(a) {
        sapply(1:ny, function(y) {
          calc_summary_F(M = M_yas[y, a, s], N = sum(N_ymars[y, 1, a, , s]),
                         CN = sum(CN_ymafrs[y, , a, , , s]), Fmax = nf * Dmodel@Fmax)
        })
      })
    })
    Z_yas <- F_yas + M_yas
    pHSP_s <- lapply(1:ns, function(s) {
      calc_HSP(yi = DCKMR@HSP_s[[s]]$yi, yj = DCKMR@HSP_s[[s]]$yj,
               N = apply(Nsp_yars[, , , s], 1:2, sum), fec = Dstock@fec_yas[, , s], Z = Z_yas[, , s])
    })
    loglike_HSP_s <- lapply(1:ns, function(s) {
      val <- 0
      if (DCKMR@HSP_s[[s]]$yi %in% y_like) {
        val <- like_CKMR(n = DCKMR@HSP_s[[s]]$n, m = DCKMR@HSP_s[[s]]$m, p = pHSP_s[[s]], type = DCKMR@CKMR_like)
      }
      return(val)
    })
  } else {
    loglike_HSP_s <- 0
  }

  ## Tag
  tag_ymarrs <- Dtag@tag_ymarrs
  if (any(tag_ymarrs > 0, na.rm = TRUE)) {
    tag_ymarrs <- OBS(tag_ymarrs)
    tagpred_ymarrs <- array(0, dim(tag_ymarrs))
    tagpred_ymarrs[] <- sapply2(1:nrow(Dtag@tag_aa), function(aa) {
      a1 <- which(Dtag@tag_aa[aa, ] > 0)[1]
      sapply2(1:nm, function(m) {
        mprev <- ifelse(m == 1, nm, m - 1)
        sapply2(1:nrow(Dtag@tag_yy), function(yy) {
          y1 <- which(Dtag@tag_yy[yy, ] > 0)[1]
          yprev <- ifelse(y1 == 1, 1, y1 - 1)
          mov_ymarrs[yprev, mprev, a1, , , ] # rrsyma
        })
      })
    }) %>%
      aperm(c(4:6, 1:3))

    loglike_tag_mov_ymars <- sapply2(1:ns, function(s) { # Likelihood of where the fish are going
      sapply2(1:nr, function(rf) {
        sapply2(1:nrow(Dtag@tag_aa), function(aa) {
          sapply(1:nm, function(m) {
            sapply(1:nrow(Dtag@tag_yy), function(yy) {
              val <- 0
              if (any(Dtag@tag_yy[yy, ] %in% y_like)) {
                pred <- tagpred_ymarrs[yy, m, aa, rf, , s]
                val <- like_comp(obs = tag_ymarrs[yy, m, aa, rf, , s],
                                 pred = CondExpGt(pred, 1e-8, pred, 1e-8), type = Dtag@tag_like,
                                 N = Dtag@tagN_ymars[yy, m, aa, rf, s], theta = Dtag@tagtheta_s[s],
                                 stdev = Dtag@tagstdev_s[s])
              }
              return(val)
            })
          })
        })
      })
    })
  } else {
    loglike_tag_mov_ymars <- 0
  }

  tag_ymars <- Dtag@tag_ymars
  if (any(tag_ymars > 0, na.rm = TRUE)) {
    tag_ymars <- OBS(tag_ymars)
    loglike_tag_dist_ymas <- sapply2(1:ns, function(s) { # Likelihood of where the fish are present
      sapply2(1:nrow(Dtag@tag_aa), function(aa) {
        avec <- Dtag@tag_aa[aa, ]
        sapply(1:nm, function(m) {
          sapply(1:nrow(Dtag@tag_yy), function(yy) {
            val <- 0
            if (any(Dtag@tag_yy[yy, ] %in% y_like)) {
              pred <- apply(N_ymars[yvec, m, avec, , s, drop = FALSE], 3, sum)
              like_comp(obs = tag_ymars[yy, m, aa, , s],
                        pred = CondExpGt(pred, 1e-8, pred, 1e-8), type = Dtag@tag_like,
                        N = Dtag@tagN_ymas[yy, m, aa, s], theta = Dtag@tagtheta_s[s],
                        stdev = Dtag@tagstdev_s[s])
            }
            return(val)
          })
        })
      })
    })
  } else {
    loglike_tag_dist_ymas <- 0
  }

  loglike <- sum(loglike_Cinit_mfr) +
    sum(loglike_Cobs_ymfr) + sum(loglike_CAA_ymfr) + sum(loglike_CAL_ymfr) +
    sum(loglike_I_ymi) + sum(loglike_IAA_ymi) + sum(loglike_IAL_ymi) +
    sum(loglike_SC_ymafr) +
    Reduce(sum, loglike_POP_s) + Reduce(sum, loglike_HSP_s) +
    sum(loglike_tag_mov_ymars) +
    sum(loglike_tag_dist_ymas)

  # Priors ----
  logprior_rdev_ys <- array(0, c(ny, ns))
  logprior_rdev_ys[] <- sapply(1:ns, function(s) {
    CondExpGt(par_rdev_ys[, s], 0, dnorm(p$log_rdev_ys[, s], 0, sdr_s[s], log = TRUE), 0)
  })
  logprior_rdev_ys[1:ny > max(y_like), ] <- 0

  logprior_initrdev_as <- sapply(1:ns, function(s) {
    CondExpGt(par_initrdev_as[, s], 0, dnorm(p$log_initrdev_as[, s], 0, sdr_s[s], log = TRUE), 0)
  })

  if (nr > 1 && "mov_g_ymars" %in% random) {
    sdg_s <- sapply2(1:ns, function(s) conv_Sigma(sigma = exp(p$log_sdg_rs[, s]), lower_diag = p$t_rhog_rs[, s]))
    par_g_ymars <- !is.na(map$mov_g_ymars) & !duplicated(map$mov_g_ymars, MARGIN = 0)

    logprior_dist_ymas <- array(0, c(ny, nm, na, ns))
    logprior_dist_ymas[y_like, , , ] <- sapply2(1:ns, function(s) {
      sapply2(1:na, function(a) {
        sapply(1:nm, function(m) {
          sapply(y_like, function(y) dmvnorm(p$mov_g_ymars[y, m, a, , s], mu = 0, Sigma = sdg_s[, , s], log = TRUE))
        })
      })
    })
  } else {
    logprior_dist_ymas <- 0
  }

  if (length(Dmodel@prior)) {
    logprior_par <- sapply(Dmodel@prior, function(x, p) eval(parse(text = x)), p = p)
  } else {
    logprior_par <- 0
  }

  logprior <- sum(logprior_par) + sum(logprior_initrdev_as) + sum(logprior_rdev_ys) + sum(logprior_dist_ymas)

  # Objective function ----
  fn <- -1 * (logprior + loglike) + penalty

  # Report out variables ----
  ## Parameters ----
  ADREPORT(R0_s)
  ADREPORT(h_s)
  REPORT(R0_s)
  REPORT(h_s)
  REPORT(kappa_s)
  REPORT(SB0_s)
  REPORT(sralpha_s)
  REPORT(srbeta_s)
  REPORT(sdr_s)
  REPORT(phi_s)

  REPORT(selconv_pf)
  if (length(Dmodel@lmid)) REPORT(sel_lf)
  REPORT(q_fs)

  if (ni > 0) {
    REPORT(selconv_pi)
    if (length(Dmodel@lmid)) REPORT(sel_li)

    ADREPORT(q_i)
  }

  REPORT(F_ymfr)

  REPORT(M_yas)
  REPORT(mat_yas)

  ## Initial (first year, first season) calculations ----
  REPORT(initNPR0_yars)
  REPORT(initRdev_as)
  if (any_Cinit) {
    REPORT(initF_mfr)
    REPORT(initZ_mars)
    REPORT(initNPR_yars)
    REPORT(initR_s)
    REPORT(initphi_s)
    REPORT(initCN_mafrs)
    REPORT(initCB_mfrs)
  }

  ## Population arrays ----
  REPORT(N_ymars)
  REPORT(S_yrs)
  REPORT(R_ys)
  REPORT(Rdev_ys)
  REPORT(F_ymars)
  REPORT(Z_ymars)
  REPORT(B_ymrs)
  if (nr > 1) {
    REPORT(mov_ymarrs)
    REPORT(recdist_rs)
  }

  ## Fishery arrays ----
  REPORT(sel_ymafs)
  REPORT(F_ymafrs)
  REPORT(CN_ymafrs)
  if (any_CAL) REPORT(CN_ymlfrs)
  REPORT(CB_ymfrs)
  REPORT(VB_ymfrs)
  if (any_SC) REPORT(SCpred_ymafrs)

  ## Index of abundance arrays ----
  if (ni > 0) {
    REPORT(sel_ymais)
    REPORT(I_ymi)
    if (any_IAA) REPORT(IN_ymais)
    if (any_IAL) REPORT(IN_ymlis)
    REPORT(q_i)
  }

  ## Tag ----
  if (any(tag_ymarrs > 0, na.rm = TRUE)) REPORT(tagpred_ymarrs)

  ## CKMR ----
  if (length(DCKMR@POP_s)) REPORT(pPOP_s)

  if (length(DCKMR@HSP_s)) {
    REPORT(F_yas)
    REPORT(Z_yas)
    REPORT(pHSP_s)
  }

  ## Objective function values ----
  REPORT(loglike)
  REPORT(logprior)
  REPORT(penalty)
  REPORT(fn)

  if (any_Cinit) REPORT(loglike_Cinit_mfr)
  if (Dmodel@condition == "F") REPORT(loglike_Cobs_ymfr)
  if (any_CAA) REPORT(loglike_CAA_ymfr)
  if (any_CAL) REPORT(loglike_CAL_ymfr)

  if (ni > 0) {
    REPORT(loglike_I_ymi)
    if (any_IAA) REPORT(loglike_IAA_ymi)
    if (any_IAL) REPORT(loglike_IAL_ymi)
  }

  if (any_SC) REPORT(loglike_SC_ymafr)

  if (length(DCKMR@POP_s)) REPORT(loglike_POP_s)
  if (length(DCKMR@HSP_s)) REPORT(loglike_HSP_s)

  if (any(tag_ymarrs > 0, na.rm = TRUE)) REPORT(loglike_tag_mov_ymars)
  if (any(tag_ymars > 0, na.rm = TRUE)) REPORT(loglike_tag_dist_ymas)

  if (length(Dmodel@prior)) REPORT(logprior_par)

  if (is.null(map$log_initrdev_as) || any(!is.na(map$log_initrdev_as))) REPORT(logprior_initrdev_as)
  if (is.null(map$log_rdev_ys) || any(!is.na(map$log_rdev_ys))) REPORT(logprior_rdev_ys)
  if (nr > 1 && "mov_g_ymars" %in% random) REPORT(logprior_dist_ymas)

  return(fn)
}

#' @importFrom methods slotNames
getAllS4 <- function (..., warn = TRUE) {
  fr <- parent.frame()
  dots <- list(...)

  for(i in 1:length(dots)) {
    nm <- slotNames(dots[[i]])
    for (j in nm) {
      if (warn && !is.null(fr[[j]])) warning("Object '", j, "' already defined")
      fr[[j]] <- slot(dots[[i]], j)
    }
  }
  invisible(NULL)
}
