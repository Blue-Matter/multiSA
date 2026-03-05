
#' @import RTMB
#' @importFrom dplyr %>%
NULL

#' @importFrom methods callNextMethod
setMethod(
  "sapply", signature(X="integer"),
  function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    ans <- callNextMethod()
    if (is.complex(ans))
      class(ans) <- "advector"
    ans
  }
)


#' `sapply2` function
#'
#' An alternate `sapply` function with argument `simplify = "array"` for convenience.
#'
#' @param X,FUN,...,USE.NAMES Same arguments as [sapply()]
#' @returns Output of `simplify2array()`, typically an array
#' @export
#' @keywords internal
sapply2 <- function(X, FUN, ..., USE.NAMES = TRUE) {
  sapply(X, FUN, ..., simplify = "array", USE.NAMES = USE.NAMES)
}

#' @name AD
#' @title Additional methods for AD types
#'
#' @description Methods for RTMB AD class
NULL

#' @describeIn AD Matrix product function implemented for mixed AD and non-AD objects with `colSums(x * y)`. See \link[RTMB]{ADmatrix}.
#'
#' @param x AD object
#' @param y Non-AD matrix
#' @aliases %*%,ad,matrix-method matmul
#' @keywords internal
setMethod("%*%",
          signature("ad", "matrix"),
          function(x, y) {
            colSums(x * y)
          })



show <- function(object) methods::show(object)
setMethod("show",
          signature = "MSAassess",
          function(object) {

            res <- list()
            res$npar <- paste("Number of parameters:", length(object@obj$par))

            cat(res$npar)
            if (length(object@SD) > 1 && !is.null(object@SD$gradient.fixed)) {
              gr <- abs(object@SD$gradient.fixed)
              gr_max <- ifelse(all(is.na(gr)), NA, round(max(gr, na.rm = TRUE), 4))

              res$max_grad <- paste("Maximum gradient:", gr_max)
            } else {
              res$max_grad <- "Run model to view gradient report"
            }
            cat("\n", res$max_grad)

            if (length(object@SD) > 1 && !is.null(object@SD$gradient.fixed)) {
              gr_na <- is.na(gr)
              if (sum(gr_na)) {
                if (sum(gr_na) < length(gr_na)) {
                  gr_names <- make_unique_names(object@SD)[gr_na]

                  res$grad_na <- paste(gr_names, collapse = ", ")

                  cat("\nParameters with gradient = NA:")
                } else {
                  res$grad_na <- "Gradient of NA for all parameters"
                }
                cat("\n", res$grad_na)
              }

              gr_large <- !is.na(gr) & gr > 0.1
              if (sum(gr_large)) {
                gr_report <- gr[gr_large]
                gr_names <- make_unique_names(object@SD)[gr_large]

                cat("\nParameters with large gradients (> 0.1):\n")

                gr_order <- order(gr_report, decreasing = TRUE)

                x <- data.frame(
                  Estimate = round(object@SD$par.fixed[gr_large], 4),
                  Gradient = round(gr_report, 4)
                )
                rownames(x) <- gr_names
                res$grad_large <- x[gr_order, ]
                print(res$grad_large)
              }
            }

            if (length(object@SD) > 1 && !is.null(object@SD$env$hessian)) {
              h <- object@SD$env$hessian

              res$det_h <- signif(det(h), 5)
              cat(paste("\nDeterminant of Hessian:"), res$det_h)

              zero_rows <- apply(h, 1, function(x) all(x == 0, na.rm = TRUE))
              na_rows <- apply(h, 1, function(x) all(is.na(x)))
              if (any(zero_rows)) {
                cat("\nParameters with all zeros in Hessian:\n")
                par_zero <- names(zero_rows)[zero_rows]
                res$hess_zero <- par_zero
                for(i in par_zero) cat(i, "\n")
              }
              if (any(na_rows)) {
                cat("\nParameters with all NAs in Hessian:\n")
                par_na <- names(na_rows)[na_rows]
                res$hess_na <- par_na
                for(i in par_na) cat(i, "\n")
              }
            }

            invisible(res)
          })
