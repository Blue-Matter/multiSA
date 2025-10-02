
#' Check dimensions and inputs in MSAdata object
#'
#' Ensures that data inputs are of proper dimension. Whenever possible, default values are added to missing items.
#'
#' @param MSAdata S4 object containing data inputs. See [MSAdata-class]
#' @param silent Logical, whether or not to report default values to the console
#' @returns An updated [MSAdata-class] object.
#' @seealso [MSAdata-class]
#' @export
check_data <- function(MSAdata, silent = FALSE) {
  MSAdata@Dmodel <- check_Dmodel(MSAdata@Dmodel, MSAdata@Dfishery@nf, silent)
  MSAdata@Dstock <- check_Dstock(MSAdata@Dstock, MSAdata@Dmodel, silent)
  MSAdata@Dfishery <- check_Dfishery(MSAdata@Dfishery, MSAdata@Dstock, MSAdata@Dmodel, silent)
  MSAdata@Dsurvey <- check_Dsurvey(MSAdata@Dsurvey, MSAdata@Dmodel, silent)

  MSAdata@DCKMR <- check_DCKMR(MSAdata@DCKMR, MSAdata@Dmodel, silent)
  MSAdata@Dtag <- check_Dtag(MSAdata@Dtag, MSAdata@Dmodel, silent)
  MSAdata@Dlabel <- check_Dlabel(MSAdata@Dlabel, MSAdata@Dmodel, MSAdata@Dfishery, MSAdata@Dsurvey, silent)
  return(MSAdata)
}


check_Dmodel <- function(Dmodel, nf, silent = FALSE) {
  ch <- as.character(substitute(Dmodel))
  if (length(ch) > 1) ch <- "Dmodel"

  if (!length(Dmodel@ny)) stop("Need ", ch, "@ny")
  if (!length(Dmodel@nm)) {
    if (!silent) message("Setting ", ch, "@nm to 1")
    Dmodel@nm <- 1L
  }
  if (!length(Dmodel@na)) stop("Need ", ch, "@na")
  if (!length(Dmodel@nl)) {
    if (!silent) message_info("No length bins in the model.")
    Dmodel@nl <- 0L
  }
  if (!length(Dmodel@nr)) {
    if (!silent) message("Creating one region model.")
    Dmodel@nr <- 1L
  }
  if (!length(Dmodel@ns)) {
    if (!silent) message("Creating one stock model.")
    Dmodel@ns <- 1L
  }
  if (Dmodel@nl > 0) {
    if (length(Dmodel@lbin) != Dmodel@nl+1) stop("Need nl + 1 vector for ", ch, "@lbin")
    if (length(Dmodel@lmid) != Dmodel@nl) stop("Need nl vector for ", ch, "@lmid")
  }
  if (!length(Dmodel@Fmax)) {
    if (!silent) message("Setting ", ch, "@Fmax to 3 (per season)")
    Dmodel@Fmax <- 3
  }
  if (!length(Dmodel@y_phi)) {
    if (!silent) message("Setting ", ch, "@y_phi to year 1")
    Dmodel@y_phi <- 1L
  }
  if (!length(Dmodel@scale_s)) {
    if (!silent) message("Setting ", ch, "@scale_s = 1 for all stocks")
    Dmodel@scale_s <- rep(1, Dmodel@ns)
  } else if (length(Dmodel@scale_s) != Dmodel@ns) {
    stop(ch, "@scale_s needs to be length ", Dmodel@ns)
  }
  if (!length(Dmodel@nyinit)) {
    if (Dmodel@nm == 1 && Dmodel@nr == 1) {
      Dmodel@nyinit <- 1
    } else {
      Dmodel@nyinit <- ceiling(1.5 * Dmodel@na)
      if (!silent) message("Setting ", ch, "@nyinit = ", Dmodel@nyinit)
    }
  }
  if (!length(Dmodel@condition)) {
    Dmodel@condition <- "F"
    if (!silent) message("Setting ", ch, "@condition = F")
  }
  if (Dmodel@condition == "catch" && !length(Dmodel@nitF)) {
    if (!silent) message("Setting ", ch, "@nitF to 5")
    Dmodel@nitF <- 5
  }
  if (Dmodel@condition == "F") {
    if (!length(Dmodel@y_Fmult_f)) {
      Dmodel@y_Fmult_f <- rep(round(0.5 * Dmodel@ny), nf)
      if (!silent) message("Setting ", ch, "@y_Fmult_f = ", paste(Dmodel@y_Fmult_f, collapse = ", "))
    }
    if (!length(Dmodel@m_Fmult_f)) {
      Dmodel@m_Fmult_f <- rep(1L, nf)
      if (!silent && Dmodel@nm > 1) message("Setting ", ch, "@m_Fmult_f = ", paste(Dmodel@m_Fmult_f, collapse = ", "))
    }
    if (!length(Dmodel@r_Fmult_f)) {
      Dmodel@r_Fmult_f <- rep(1L, nf)
      if (!silent && Dmodel@nr > 1) message("Setting ", ch, "@r_Fmult_f = ", paste(Dmodel@r_Fmult_f, collapse = ", "))
    }
  }

  return(Dmodel)
}

check_Dstock <- function(Dstock, Dmodel, silent = FALSE) {
  ny <- Dmodel@ny
  nm <- Dmodel@nm
  na <- Dmodel@na
  nl <- Dmodel@nl
  nr <- Dmodel@nr
  ns <- Dmodel@ns

  ch <- as.character(substitute(Dstock))
  if (length(ch) > 1) ch <- "Dstock"

  if (!length(Dstock@m_spawn) || nm == 1L) {
    if (!silent && nm > 1) message("Setting ", ch, "@m_spawn = 1")
    Dstock@m_spawn <- 1L
  } else if (Dstock@m_spawn > nm) {
    stop(ch, "@m_spawn cannot be greater than Dmodel@nm")
  }
  if (!length(Dstock@m_rec) || nm == 1L) {
    if (!silent && Dmodel@nm > 1) message("Setting ", ch, "@m_rec = 1")
    Dstock@m_rec <- 1L
  } else if (Dstock@m_rec > nm) {
    stop(ch, "@m_rec cannot be greater than Dmodel@nm")
  }

  if (nl > 0) {
    if (length(Dstock@LAK_ymals)) {
      dim_LAK <- dim(Dstock@LAK_ymals) == c(ny, nm, na, nl, ns)
      if (!all(dim_LAK)) stop("dim(", ch, "@LAK_ymals) needs to be: ", c(ny, nm, na, nl, ns) %>% paste(collapse = ", "))
    } else {

      dim_len <- dim(Dstock@len_ymas) == c(ny, nm, na, ns)
      if (!all(dim_len)) stop("dim(", ch, "@len_ymas) needs to be: ", c(ny, nm, na, ns) %>% paste(collapse = ", "))

      dim_sdlen <- dim(Dstock@sdlen_ymas) == c(ny, nm, na, ns)
      if (!all(dim_sdlen)) stop("dim(", ch, "@sdlen_ymas) needs to be: ", c(ny, nm, na, ns) %>% paste(collapse = ", "))

      if (!silent) message("Calculating ", ch, "@LAK_ymals array")
      Dstock@LAK_ymals <- sapply2(1:ns, function(s) {
        sapply2(1:nm, function(m) {
          sapply2(1:ny, function(y) calc_LAK(Dstock@len_ymas[y, m, , s], Dstock@sdlen_ymas[y, m, , s], Dmodel@lbin))
        })
      }) %>% aperm(c(3, 4, 1, 2, 5))
    }
  }

  if (!length(Dstock@matd_yas)) {
    dim_mat <- dim(Dstock@matd_yas) == c(ny, na, ns)
    if (!all(dim_mat)) stop("dim(", ch, "@matd_yas) needs to be: ", c(ny, na, ns) %>% paste(collapse = ", "))
  }

  dim_swt <- dim(Dstock@swt_ymas) == c(ny, nm, na, ns)
  if (!all(dim_swt)) stop("dim(", ch, "@swt_ymas) needs to be: ", c(ny, nm, na, ns) %>% paste(collapse = ", "))

  if (!length(Dstock@fec_yas)) {
    if (!silent) message("Setting fecundity to stock weight at age at season of spawning")
    Dstock@fec_yas <- array(Dstock@swt_ymas[, Dstock@m_spawn, , ], c(ny, na, ns))
  } else {
    dim_fec <- dim(Dstock@fec_yas) == c(ny, na, ns)
    if (!all(dim_fec)) stop("dim(", ch, "fec_yas) needs to be: ", c(ny, na, ns) %>% paste(collapse = ", "))
  }

  if (!length(Dstock@Md_yas)) {
    dim_M <- dim(Dstock@Md_yas) == c(ny, na, ns)
    if (!all(dim_M)) stop("dim(", ch, "@Md_yas) needs to be: ", c(ny, na, ns) %>% paste(collapse = ", "))
  }

  if (length(Dstock@SRR_s) != ns) stop("SRR_s needs to be length ", ns)
  if (!length(Dstock@delta_s)) {
    Dstock@delta_s <- rep(0, ns)
  } else if (length(Dstock@delta_s) != ns) {
    stop(ch, "@delta_s needs to be length ", ns)
  }
  if (!length(Dstock@presence_rs)) {
    Dstock@presence_rs <- matrix(TRUE, nr, ns)
  } else if (any(dim(Dstock@presence_rs) != c(nr, ns))) {
    stop("dim(", ch, "@presence_rs) needs to be: ", c(nr, ns) %>% paste(collapse = ", "))
  }
  if (!length(Dstock@natal_rs)) {
    Dstock@natal_rs <- matrix(1, nr, ns)
  } else if (any(dim(Dstock@natal_rs) != c(nr, ns))) {
    stop("dim(", ch, "@natal_rs) needs to be: ", c(nr, ns) %>% paste(collapse = ", "))
  }

  return(Dstock)
}

check_Dfishery <- function(Dfishery, Dstock, Dmodel, silent = FALSE) {
  ny <- Dmodel@ny
  nm <- Dmodel@nm
  na <- Dmodel@na
  nl <- Dmodel@nl
  nr <- Dmodel@nr
  ns <- Dmodel@ns

  ch <- as.character(substitute(Dfishery))
  if (length(ch) > 1) ch <- "Dfishery"

  if (!length(Dfishery@nf)) stop("Need ", ch, "@nf")
  nf <- Dfishery@nf

  dim_Cobs <- dim(Dfishery@Cobs_ymfr) == c(ny, nm, nf, nr)
  if (!all(dim_Cobs)) stop("dim(", ch, "@Cobs_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))

  if (Dmodel@condition == "F") {
    if (!length(Dfishery@Csd_ymfr)) {
      Dfishery@Csd_ymfr <- array(0.01, c(ny, nm, nf, nr))
    } else {
      dim_Csd <- dim(Dfishery@Csd_ymfr) == c(ny, nm, nf, nr)
      if (!all(dim_Csd)) stop("dim(", ch, "@Csd_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))
    }
  }

  if (!length(Dfishery@fwt_ymafs)) {
    if (!silent) message("Setting fishery weight at age to stock weight at age")
    Dfishery@fwt_ymafs <- sapply2(1:ns, function(s) {
      sapply2(1:nf, function(f) array(Dstock@swt_ymas[, , , s], c(ny, nm, na)))
    })
  } else {
    dim_fwt <- dim(Dfishery@fwt_ymafs) == c(ny, nm, na, nf, ns)
    if (!all(dim_fwt)) stop("dim(", ch, "@fwt_ymafs) needs to be: ", c(ny, nm, na, nf, ns) %>% paste(collapse = ", "))
  }

  if (length(Dfishery@CAAobs_ymafr) || length(Dfishery@CALobs_ymlfr)) {
    if (length(Dfishery@fcomp_like)) {
      Dfishery@fcomp_like <- match.arg(Dfishery@fcomp_like, choices = eval(formals(like_comp)$type))
    } else {
      if (!silent) message("Setting ", ch, "@fcomp_like = \"multinomial\"")
      Dfishery@fcomp_like <- "multinomial"
    }
  }

  if (length(Dfishery@CAAobs_ymafr)) {
    dim_CAA <- dim(Dfishery@CAAobs_ymafr) == c(ny, nm, na, nf, nr)
    if (!all(dim_CAA)) stop("dim(", ch, "@CAAobs_ymafr) needs to be: ", c(ny, nm, na, nf, nr) %>% paste(collapse = ", "))

    if (!length(Dfishery@CAAN_ymfr)) {
      if (Dfishery@fcomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@CAAN_ymfr from ", ch, "@CAAobs_ymafr")
      }
      Dfishery@CAAN_ymfr <- apply(Dfishery@CAAobs_ymafr, c(1, 2, 4, 5), sum)
    } else {
      dim_CAAN <- dim(Dfishery@CAAN_ymfr) == c(ny, nm, nf, nr)
      if (!all(dim_CAAN)) stop("dim(", ch, "@CAAN_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))
    }

    if (!length(Dfishery@CAAtheta_f)) {
      if (grepl("ddirmult", Dfishery@fcomp_like) && !silent) message("Setting ", ch, "@CAAtheta_f to 1 for all fleets")
      Dfishery@CAAtheta_f <- rep(1, nf)
    } else if (length(Dfishery@CAAtheta_f) == 1) {
      Dfishery@CAAtheta_f <- rep(Dfishery@CAAtheta_f, nf)
    } else if (length(Dfishery@CAAtheta_f) != nf) {
      stop("Vector ", ch, "@CAAtheta_f needs to be length ", nf)
    }
  }

  if (length(Dfishery@CALobs_ymlfr)) {
    if (!nl) stop("Fishery length composition detected but no length bins found in Dmodel@nl")

    dim_CAL <- dim(Dfishery@CALobs_ymlfr) == c(ny, nm, nl, nf, nr)
    if (!all(dim_CAL)) stop("dim(", "@CALobs_ymlfr) needs to be: ", c(ny, nm, nl, nf, nr) %>% paste(collapse = ", "))

    if (!length(Dfishery@CALN_ymfr)) {
      if (Dfishery@fcomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@CALN_ymfr from ", ch, "@CALobs_ymlfr")
      }
      Dfishery@CALN_ymfr <- apply(Dfishery@CALobs_ymlfr, c(1, 2, 4, 5), sum)
    } else {
      dim_CALN <- dim(Dfishery@CALN_ymfr) == c(ny, nm, nf, nr)
      if (!all(dim_CALN)) stop("dim(", ch, "@CALN_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))
    }

    if (!length(Dfishery@CALtheta_f)) {
      if (grepl("ddirmult", Dfishery@fcomp_like) && !silent) message("Setting ", ch, "@CALtheta_f to 1 for all fleets")
      Dfishery@CALtheta_f <- rep(1, nf)
    } else if (length(Dfishery@CALtheta_f) == 1) {
      Dfishery@CALtheta_f <- rep(Dfishery@CALtheta_f, nf)
    } else if (length(Dfishery@CALtheta_f) != nf) {
      stop("Vector ", ch, "@CALtheta_f needs to be length ", nf)
    }
  }

  if (!length(Dfishery@sel_block_yf)) {
    Dfishery@sel_block_yf <- matrix(1:nf, ny, nf, byrow = TRUE)
  } else {
    dim_sb <- dim(Dfishery@sel_block_yf) == c(ny, nf)
    if (!all(dim_sb)) stop("dim(", ch, "@sel_block_yf) needs to be: ", c(ny, nf) %>% paste(collapse = ", "))
  }
  nb <- max(Dfishery@sel_block_yf)
  if (length(Dfishery@sel_f) != nb) stop("Vector sel_f should be length ", nf)

  if (!length(Dfishery@Cinit_mfr)) {
    Dfishery@Cinit_mfr <- array(0, c(nm, nf, nr))
  } else {
    dim_Cinit <- dim(Dfishery@Cinit_mfr) == c(nm, nf, nr)
    if (!all(dim_Cinit)) stop("dim(", ch, "@Cinit_mfr) needs to be: ", c(nm, nf, nr) %>% paste(collapse = ", "))
  }

  if (length(Dfishery@SC_ymafrs)) {
    dim_SC <- dim(Dfishery@SC_ymafrs)
    if (length(dim_SC) != 6) stop(ch, "@SC_ymafrs should be a six dimensional array")

    if (any(dim_SC[c(1, 2, 5, 6)] != c(ny, nm, nr, ns))) {
      stop("dim(", ch, "@SC_ymafrs) should be ", c(ny, nm, dim_SC[3:4], nr, ns) %>% paste(collapse = ", "))
    }

    if (dim_SC[3] == na) {
      if (!length(Dfishery@SC_aa)) Dfishery@SC_aa <- diag(1, na)
    } else if (any(dim(Dfishery@SC_aa) != c(dim_SC[3], na))) {
      stop("dim(Dfishery@SC_aa) should be: ", c(Dfishery@dim_SC[3], na) %>% paste(collapse = ", "))
    }

    if (dim_SC[4] == nf) {
      if (!length(Dfishery@SC_ff)) Dfishery@SC_ff <- diag(1, nf)
    } else if (any(dim(Dfishery@SC_ff) != c(dim_SC[4], nf))) {
      stop("dim(", ch, "@SC_ff) should be: ", c(dim_SC[4], nf) %>% paste(collapse = ", "))
    }

    if (!length(Dfishery@SC_like)) {
      if (!silent) message("Setting ", ch, "@SC_like = \"multinomial\"")
      Dfishery@SC_like <- "multinomial"
    } else {
      Dfishery@SC_like <- match.arg(Dfishery@SC_like, choices = eval(formals(like_comp)$type))
    }

    if (!length(Dfishery@SCN_ymafr)) {
      if (Dfishery@SC_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@SCN_ymafr from ", ch, "@SC_ymafrs")
      }
      Dfishery@SCN_ymafr <- apply(Dfishery@SC_ymafrs, 1:5, sum)
    } else {
      dim_SCN <- dim(Dfishery@SCN_ymafr) == c(ny, nm, dim_SC[3], dim_SC[4], nr)
      if (!all(dim_SCN)) stop("dim(", ch, "@SCN_ymafr) needs to be: ", c(ny, nm, dim_SC[3], dim_SC[4], nr) %>% paste(collapse = ", "))
    }

    if (!length(Dfishery@SCtheta_f)) {
      if (grepl("ddirmult", Dfishery@SC_like) && !silent) message("Setting ", ch, "@SCtheta_f to 1 for all fleets")
      Dfishery@SCtheta_f <- rep(1, dim_SC[4])
    } else if (length(Dfishery@SCtheta_f) == 1) {
      Dfishery@SCtheta_f <- rep(Dfishery@SCtheta_f, dim_SC[4])
    } else if (length(Dfishery@SCtheta_f) != dim_SC[4]) {
      stop("Vector ", ch, "@SCtheta_f needs to be length ", dim_SC[4])
    }

    if (!length(Dfishery@SCstdev_ymafrs)) {
      if (grepl("log", Dfishery@SC_like) && !silent) message("Setting ", ch, "@SCstdev_ymafrs to 0.1 for all fleets")
      Dfishery@SCstdev_ymafrs <- array(0.1, c(ny, nm, dim_SC[3], dim_SC[4], nr, ns))
    } else if (length(Dfishery@SCstdev_ymafrs) == 1) {
      Dfishery@SCstdev_ymafrs <- array(Dfishery@SCstdev_ymafrs, c(ny, nm, dim_SC[3], dim_SC[4], nr, ns))
    } else if (any(dim(Dfishery@SCstdev_ymafrs) != c(ny, nm, dim_SC[3], dim_SC[4], nr, ns))) {
      stop("dim(", ch, "@SCstdev_ymafrs) needs to be: ", c(ny, nm, dim_SC[3], dim_SC[4], nr, ns) %>% paste(collapse = ", "))
    }
  }
  return(Dfishery)
}

check_Dsurvey <- function(Dsurvey, Dmodel, silent = FALSE) {
  ny <- Dmodel@ny
  nm <- Dmodel@nm
  na <- Dmodel@na
  nl <- Dmodel@nl
  nr <- Dmodel@nr
  ns <- Dmodel@ns

  ch <- as.character(substitute(Dsurvey))
  if (length(ch) > 1) ch <- "Dsurvey"

  if (!length(Dsurvey@ni)) {
    if (length(Dsurvey@Iobs_ymi) || length(Dsurvey@IAAobs_ymai) || length(Dsurvey@IALobs_ymli)) {
      stop("Need ", ch, "@ni")
    } else {
      if (!silent) message("Setting ", ch, "@ni to zero")
      Dsurvey@ni <- 0
    }
  }
  ni <- Dsurvey@ni

  if (ni > 0) {
    dim_Iobs <- dim(Dsurvey@Iobs_ymi) == c(ny, nm, ni)
    if (!all(dim_Iobs)) stop("dim(", ch, "@Iobs_ymi) needs to be: ", c(ny, nm, ni) %>% paste(collapse = ", "))

    if (!length(Dsurvey@unit_i)) {
      if (!silent) message("Setting unit_i to biomass for all indices")
      Dsurvey@unit_i <- rep("B", ni)
    }

    if (length(Dsurvey@IAAobs_ymai) || length(Dsurvey@IALobs_ymli)) {
      if (length(Dsurvey@icomp_like)) {
        Dsurvey@icomp_like <- match.arg(Dsurvey@icomp_like, choices = eval(formals(like_comp)$type))
      } else {
        if (!silent) message("Setting ", ch, "@icomp_like = \"multinomial\"")
        Dsurvey@icomp_like <- "multinomial"
      }
    }

    if (length(Dsurvey@IAAobs_ymai)) {
      dim_IAA <- dim(Dsurvey@IAAobs_ymai) == c(ny, nm, na, ni)
      if (!all(dim_IAA)) stop("dim(IAAobs_ymai) needs to be: ", c(ny, nm, na, ni) %>% paste(collapse = ", "))

      if (!length(Dsurvey@IAAN_ymi)) {
        if (Dsurvey@icomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
          message("Setting ", ch, "@IAAN_ymi from ", ch, "@IAAobs_ymai")
        }
        Dsurvey@IAAN_ymi <- apply(Dsurvey@IAAobs_ymai, c(1, 2, 4), sum)
      } else {
        dim_IAAN <- dim(Dsurvey@IAAN_ymi) == c(ny, nm, ni)
        if (!all(dim_IAAN)) stop("dim(", ch, "@IAAN_ymi) needs to be: ", c(ny, nm, ni) %>% paste(collapse = ", "))
      }

      if (!length(Dsurvey@IAAtheta_i)) {
        if (grepl("ddirmult", Dsurvey@icomp_like) && !silent) message("Setting ", ch, "@IAAtheta_i to 1 for all indices")
        Dsurvey@IAAtheta_i <- rep(1, ni)
      } else if (length(Dsurvey@IAAtheta_i) == 1) {
        Dsurvey@IAAtheta_i <- rep(Dsurvey@IAAtheta_i, ni)
      } else if (length(Dsurvey@IAAtheta_i) != ni) {
        stop("Vector ", ch, "@IAAtheta_i needs to be length ", ni)
      }
    }

    if (length(Dsurvey@IALobs_ymli)) {
      if (!nl) stop("Index length composition detected but no length bins found in Dmodel@nl")

      dim_IAL <- dim(Dsurvey@IALobs_ymli) == c(ny, nm, nl, ni)
      if (!all(dim_IAL)) stop("dim(", ch, "@IALobs_ymli) needs to be: ", c(ny, nm, nl, ni) %>% paste(collapse = ", "))

      if (!length(Dsurvey@IALN_ymi)) {
        if (Dsurvey@icomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
          message("Setting ", ch, "@IALN_ymi from ", ch, "@IALobs_ymli")
        }
        Dsurvey@IALN_ymi <- apply(Dsurvey@IALobs_ymli, c(1, 2, 4), sum)
      } else {
        dim_IALN <- dim(Dsurvey@IALN_ymi) == c(ny, nm, ni)
        if (!all(dim_IALN)) stop("dim(", ch, "@IALN_ymi) needs to be: ", c(ny, nm, ni) %>% paste(collapse = ", "))
      }

      if (!length(Dsurvey@IALtheta_i)) {
        if (grepl("ddirmult", Dsurvey@icomp_like) && !silent) message("Setting ", ch, "@IALtheta_i to 1 for all indices")
        Dsurvey@IALtheta_i <- rep(1, ni)
      } else if (length(Dsurvey@IALtheta_i) == 1) {
        Dsurvey@IALtheta_i <- rep(Dsurvey@IALtheta_i, ni)
      } else if (length(Dsurvey@IALtheta_i) != ni) {
        stop("Vector ", ch, "@IALtheta_i needs to be length ", ni)
      }
    }

    if (!length(Dsurvey@samp_irs)) {
      if (nr > 1 || ns > 1) {
        if (!silent) message("Setting ", ch, "@samp_irs = 1. All indices operate in all regions and sample all stocks")
      }
      Dsurvey@samp_irs <- array(1, c(ni, nr, ns))
    } else {
      dim_samp <- dim(Dsurvey@samp_irs) == c(ni, nr, ns)
      if (!all(dim_samp)) stop("dim(", ch, "@samp_irs) needs to be: ", c(ni, nr, ns) %>% paste(collapse = ", "))
    }

    if (length(Dsurvey@sel_i) == 1) {
      Dsurvey@sel_i <- rep(Dsurvey@sel_i, ni)
    } else if(length(Dsurvey@sel_i) != ni) {
      stop("Vector ", ch, "@sel_i needs to be length ", ni)
    }

    if (!length(Dsurvey@delta_i)) {
      if (!silent) message("Setting ", ch, "@delta_i = 0 for all indices (survey timing within time step)")
      Dsurvey@delta_i <- rep(0, ni)
    } else if (length(Dsurvey@delta_i) == 1) {
      Dsurvey@delta_i <- rep(Dsurvey@delta_i, ni)
    } else if(length(Dsurvey@delta_i) != ni) {
      stop("Vector ", ch, "@delta_i needs to be length ", ni)
    }
  }

  return(Dsurvey)
}

check_DCKMR <- function(DCKMR, Dmodel, silent = FALSE) {
  ns <- Dmodel@ns

  if (length(DCKMR@POP_s)) {
    vars_POP <- c("a", "t", "y", "n", "m")
    check_POP <- sapply(1:ns, function(s) {
      if (nrow(DCKMR@POP_s[[s]])) {
        all(vars_POP %in% names(DCKMR@POP_s[[s]]))
      } else {
        TRUE
      }
    })
    if (any(!check_POP)) {
      stop("Missing columns in close-kin POP data frames. See: help(\"MSAdata-class\")")
    }
  }

  if (length(DCKMR@HSP_s)) {
    vars_HSP <- c("yi", "yj", "n", "m")
    check_HSP <- sapply(1:ns, function(s) {
      if (nrow(DCKMR@HSP_s[[s]])) {
        all(vars_HSP %in% names(DCKMR@HSP_s[[s]]))
      } else {
        TRUE
      }
    })
    if (any(!check_HSP)) {
      stop("Missing columns in close-kin HSP data frames. See: help(\"MSAdata-class\")")
    }
  }

  if (length(DCKMR@POP_s) || length(DCKMR@HSP_s)) {
    if (!length(DCKMR@CKMR_like)) {
      if (!silent) message("Setting close-kin likelihood to \"binomial\"")
      DCKMR@CKMR_like <- "binomial"
    }
  }
  return(DCKMR)
}

check_Dtag <- function(Dtag, Dmodel, silent = FALSE) {
  ny <- Dmodel@ny
  nm <- Dmodel@nm
  na <- Dmodel@na
  nr <- Dmodel@nr
  ns <- Dmodel@ns

  ch <- as.character(substitute(Dtag))
  if (length(ch) > 1) ch <- "Dtag"

  if (length(Dtag@tag_ymarrs) && length(Dtag@tag_ymars)) {
    stop("Tag data found in both Dtag@tag_ymarrs (tag movement) and Dtag@ymars (tag distribution). Only one tag setup can be used at this time.")
  }

  if (length(Dtag@tag_ymarrs) || length(Dtag@tag_ymars)) {
    if (length(Dtag@tag_like)) {
      Dtag@tag_like <- match.arg(Dtag@tag_like, choices = eval(formals(like_comp)$type))
    } else {
      if (!silent) message("Setting ", ch, "@tag_like = \"multinomial\"")
      Dtag@tag_like <- "multinomial"
    }

    if (!length(Dtag@tagtheta_s)) {
      if (grepl("ddirmult", Dtag@tag_like) && !silent) message("Setting ", ch, "@tagtheta_s to 1 for all stocks")
      Dtag@tagtheta_s <- rep(1, ns)
    } else if (length(Dtag@tagtheta_s) == 1) {
      Dtag@tagtheta_s <- rep(Dtag@tagtheta_s, ns)
    } else if (length(Dtag@tagtheta_s) != ns) {
      stop("Vector ", ch, "@tagtheta_s needs to be length ", ns)
    }

    if (!length(Dtag@tagstdev_s)) {
      if (grepl("log", Dtag@tag_like) && !silent) message("Setting ", ch, "@tagstdev_s to 0.1 for all fleets")
      Dtag@tagstdev_s <- rep(0.1, ns)
    } else if (length(Dtag@tagstdev_s) == 1) {
      Dtag@tagstdev_s <- rep(Dtag@tagstdev_s, ns)
    } else if (length(Dtag@tagstdev_s) != ns) {
      stop("Vector ", ch, "@tagstdev_s needs to be length ", ns)
    }
  } else { # Sets up movement estimation
    Dtag@tag_yy <- matrix(1:ny, 1, ny)
  }

  if (length(Dtag@tag_ymarrs)) {

    dim_tag1 <- dim(Dtag@tag_ymarrs)
    if (length(dim_tag1) != 6) stop("tag_ymarrs should be a six dimensional array")

    if (any(dim_tag1[c(2, 4, 5, 6)] != c(nm, nr, nr, ns))) {
      stop("dim(", ch, "@tag_ymarrs) should be ", c(dim_tag1[1], nm, dim_tag1[2], nr, nr, ns) %>% paste(collapse = ", "))
    }

    if (dim_tag1[1] == ny) {
      if (!length(Dtag@tag_yy)) Dtag@tag_yy <- diag(1, ny)
    } else if (any(dim(Dtag@tag_yy) != c(dim_tag1[1], ny))) {
      stop("dim(", ch, "@tag_yy) should be: ", c(dim_tag1[1], ny) %>% paste(collapse = ", "))
    }

    if (dim_tag1[3] == na) {
      if (!length(Dtag@tag_aa)) Dtag@tag_aa <- diag(1, na)
    } else if (any(dim(Dtag@tag_aa) != c(dim_tag1[3], na))) {
      stop("dim(", ch, "@tag_aa) should be: ", c(dim_tag1[3], na) %>% paste(collapse = ", "))
    }

    if (!length(Dtag@tagN_ymars)) {
      if (Dtag@tag_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@tagN_ymars from ", ch, "@tag_ymarrs")
      }
      Dtag@tagN_ymars <- apply(Dtag@tag_ymarrs, c(1:4, 6), sum)
    } else {
      dim_tagN <- dim(Dtag@tagN_ymars) == c(dim_tag1[1], nm, dim_tag1[3], nr, ns)
      if (!all(dim_tagN)) stop("dim(", ch, "@tagN_ymars) needs to be: ", c(dim_tag1[1], nm, dim_tag1[3], nr, ns) %>% paste(collapse = ", "))
    }

  }

  if (length(Dtag@tag_ymars)) {

    dim_tag2 <- dim(Dtag@tag_ymars)
    if (length(dim_tag2) != 5) stop("tag_ymars should be a five dimensional array")

    if (any(dim_tag2[c(1, 2, 4, 5)] != c(ny, nm, nr, ns))) {
      stop("dim(tag_ymars) should be ", c(ny, nm, dim_tag2[3], nr, ns) %>% paste(collapse = ", "))
    }

    if (dim_tag2[3] == na) {
      if (!length(Dtag@tag_aa)) Dtag@tag_aa <- diag(1, na)
    } else if (any(dim(Dtag@tag_aa) != c(dim_tag2[3], na))) {
      stop("dim(tag_aa) should be: ", c(dim_tag2[3], na) %>% paste(collapse = ", "))
    }

    if (!length(Dtag@tagN_ymas)) {
      if (Dtag@tag_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@tagN_ymas from ", ch, "@tag_ymars")
      }
      Dtag@tagN_ymas <- apply(Dtag@tag_ymars, c(1:3, 5), sum)
    } else {
      dim_tagN <- dim(Dtag@tagN_ymas) == c(dim_tag2[1], nm, dim_tag2[3], ns)
      if (!all(dim_tagN)) stop("dim(tagN_ymas) needs to be: ", c(dim_tag2[1], nm, dim_tag2[3], ns) %>% paste(collapse = ", "))
    }

  }

  return(Dtag)
}

check_Dlabel <- function(Dlabel, Dmodel, Dfishery, Dsurvey, silent = FALSE) {
  ny <- Dmodel@ny
  nm <- Dmodel@nm
  na <- Dmodel@na
  nr <- Dmodel@nr
  ns <- Dmodel@ns

  nf <- Dfishery@nf
  ni <- Dsurvey@ni

  ch <- as.character(substitute(Dlabel))
  if (length(ch) > 1) ch <- "Dlabel"

  if (!length(Dlabel@year)) {
    Dlabel@year <- 1:ny
  } else if (length(Dlabel@year) != ny) {
    stop("length(", ch, "@year) needs to be ", ny)
  }
  if (nm > 1) {
    if (!length(Dlabel@season) && nm > 1) {
      Dlabel@season <- paste("Season", 1:nm)
    } else if (length(Dlabel@season) != nm) {
      stop("length(", ch, "@season) needs to be ", nm)
    }
  }
  if (!length(Dlabel@age)) {
    Dlabel@age <- 1:na
  } else if (length(Dlabel@age) != na) {
    stop("length(", ch, "@age) needs to be ", na)
  }
  if (nr > 1) {
    if (!length(Dlabel@region) && nr > 1) {
      Dlabel@region <- paste("Region", 1:nr)
    } else if (length(Dlabel@region) != nr) {
      stop("length(", ch, "@region) needs to be ", nr)
    }
  }
  if (ns > 1) {
    if (!length(Dlabel@stock) && ns > 1) {
      Dlabel@stock <- paste("Stock", 1:ns)
    } else if (length(Dlabel@stock) != ns) {
      stop("length(", ch, "@stock) needs to be ", ns)
    }
  }
  if (!length(Dlabel@fleet)) {
    Dlabel@fleet <- paste("Fleet", 1:nf)
  } else if (length(Dlabel@fleet) != nf) {
    stop("length(", ch, "@fleet) needs to be ", nf)
  }
  if (Dsurvey@ni > 0) {
    if (!length(Dlabel@index)) {
      Dlabel@index <- paste("Index", 1:ni)
    } else if (length(Dlabel@index) != ni) {
      stop("length(", ch, "@index) needs to be ", ni)
    }
  }

  return(Dlabel)
}
