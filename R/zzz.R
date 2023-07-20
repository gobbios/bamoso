#' grooming matrices and observation effort for four macaque species
#'
#' @format a list with four items. Each item is a list itself with one grooming
#'         matrix and a vector or matrix for observation effort.
#' @source
#' \itemize{
#'  \item{\emph{Macaca assamensis}:} {Oliver Schülke and Julia Ostner}
#'  \item{\emph{M. fuscata}:} {Andrew MacIntosh}
#'  \item{\emph{M. nigra}:} {Julie Duboscq}
#'  \item{\emph{M. sylvanus}:} {Julia Fischer}
#' }
"grooming"

#' socprog example: raw association data and metrics calculated within SOCPROG.
#'
#' @references
#' \insertRef{cairns1987}{basr}
#'
#' \insertRef{whitehead2008}{basr}
#'
#' \insertRef{whitehead2015}{basr}
#'
#' @format a list with two items. First item are the raw associations,
#'         second are the different metrics returned by SOCPROG.
"socprogexample"

#' data for workflow example
#'
#'
#' @format a list with named items.
#' @examples
#' set.seed(755907)
#' workflow_example <- generate_data(n_ids = 25,
#'                                   n_beh = 1,
#'                                   behav_types = "count",
#'                                   indi_sd = 0.9,
#'                                   dyad_sd = 2.3,
#'                                   beh_intercepts = -1.5)
#' indivals <- workflow_example$input_data$indi_soc_vals
#' workflow_example$response <- rpois(n = workflow_example$standat$n_ids,
#'                                    lambda = exp(-0.5 + 1.2 * indivals))
#'
"workflow_example"
