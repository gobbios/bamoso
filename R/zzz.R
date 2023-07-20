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

#' data for kangaroo reproductive success
#' @details
#' \code{kangaroos_rawdata} contains a subset 47 individuals of the raw data
#' taken from the primary sources (see references). It is formatted as a list.
#' The first item contains as a list the raw association data (one matrix
#' per year). The second item lists the fitness outcomes and values for
#' control variables for each individual in each year.
#'
#' \code{kangaroos_subset} contains a list derived from the association data
#' and the fitness outcome taken from \code{kangaroos_rawdata}. It also contains
#' the values of the variables included in the model that were not directly
#' quantified from the association data (parity status, residual group size,
#' number of sightings). It is formatted such that it can be passed to Stan.
#'
#' @aliases kangaroos_rawdata
#'
#' @references
#' \insertRef{menz2020}{basr}
#'
#' @source
#' \insertRef{menz2020data}{basr}
#'
#' @format a list with named items to be passed to Stan.
"kangaroos_subset"
