#' @description S4 class that organizes the various data inputs for the MSA model. `MSAdata` simply inherits the slots from 6 component classes:
#' `Dmodel`, `Dstock`, `Dfishery`, `Dsurvey` `DCKMR`, and `Dtag`, where the `D`- prefix denotes an object for data inputs (or model configuration).
#'
#' @details
#' For convenience, most arrays and matrices have the associated dimensions in the variable name. For example, `Cobs_ymfr` represents
#' observed catch with the dimension following the underscore, following this template:
#'
#' \tabular{ll}{
#' `y` \tab Year \cr
#' `m` \tab Season \cr
#' `a` \tab Age \cr
#' `r` \tab Region \cr
#' `f` \tab Fishery \cr
#' `i` \tab Index \cr
#' `s` \tab Stock
#' }
#' @seealso \link{MSAdata-class} [check_data()] \link{Dmodel-class} \link{Dstock-class} \link{Dfishery-class} \link{Dsurvey-class} \link{DCKMR-class} \link{Dtag-class}
