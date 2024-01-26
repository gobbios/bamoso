#' fit dyadic relationship model to one or more interaction matrices
#'
#' @param standat stan data list (typically the result of a call
#'                 to \code{\link{make_stan_data_from_matrices}} or
#'                 \code{\link{make_stan_data_from_association}})
#' @param ... further arguments for \code{\link[cmdstanr]{sample}} (typically
#'            \code{parallel_chains}, \code{refresh},
#'            \code{iter_warmup}, \code{iter_sampling},
#'            \code{seed}, \code{adapt_delta},
#'            and/or \code{step_size})
#'
#' @return a list (of class \code{dyadicmodel}) which contains the standata
#'         list as first item and the CmdStanFit as second item
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("grooming")
#' mat <- grooming$ass$groom
#' obseff <- grooming$ass$obseff
#' standat <- make_stan_data_from_matrices(mats = list(groom = mat),
#'                                         behav_types = "count",
#'                                         obseff = list(obseff))
#' res <- sociality_model(standat = standat, parallel_chains = 4,
#'                        adapt_delta = 0.9, seed = 123)
#' summary(res)
#' }
#'
#' \dontrun{
#' x <- generate_data(n_ids = 16, n_beh = 1, indi_sd = 1,
#'                    dyad_sd = 1, beh_intercepts = 1)
#' standat <- x$standat
#' res <- sociality_model(standat = standat, parallel_chains = 4)
#' }

sociality_model <- function(standat,
                            ...) {

  modeltype <- "simple"
  if (!is.null(standat$indi_cat_pred)) {
    modeltype <- "indi_cat"
  }
  mod <- get_model(type = modeltype)
  res <- mod$sample(data = standat, ...)
  out <- list(standat = standat, mod_res = res)
  class(out) <- "dyadicmodel"
  out
}
