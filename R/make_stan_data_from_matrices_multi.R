#' prepare data for Stan for multi-group model
#'
#' @param bmats list with named(!) interaction matrices (undirected,
#'          upper triangle only)
#' @param btypes character: either \code{"count"} or \code{"prop"}
#' @param omats optional list of named matrices with observation effort
#' @param forsimeval internal, for simulations only
#' @param indidata internal, for simulations only
#' @param dyaddata internal, for simulations only
#'
#' @returns a list ready to be passed to Stan
#'
#' @examples
#' data("grooming")
#' bmats <- list(nig = grooming$nig$groom, ass = grooming$ass$groom)
#' omats <- list(nig = grooming$nig$obseff, ass = grooming$ass$obseff)
#' bamoso:::make_stan_data_from_matrices_multi(bmats, "count", omats)

make_stan_data_from_matrices_multi <- function(bmats,
                                               btypes = c("count", "prop"),
                                               omats = NULL,
                                               forsimeval = FALSE,
                                               indidata = NULL,
                                               dyaddata = NULL
                                               ) {

  if (btypes == "prop" && is.null(omats)) {
    stop("require observation effort for discrete proportion data")
  }

  if (btypes == "count" && is.null(omats)) {
    omats <- lapply(bmats, \(x){
      x <- x - x
      x[upper.tri(x)] <- 1
      x
    })
  }

  btyp_num <- NULL
  if (btypes == "count") btyp_num <- c(count = 1)
  if (btypes == "prop") btyp_num <- c(prop = 2)

  groups <- lapply(bmats, colnames)
  aux <- make_indices(grouplist = groups)
  indilist <- aux$indices_indi
  indices <- aux$indices_dyad
  indices$b <- NA
  indices$obseff <- NA

  for (i in seq_len(nrow(indices))) {
    b <- bmats[[indices$period_char[i]]]
    o <- omats[[indices$period_char[i]]]
    indices$b[i] <- b[indices$within_id1_char[i],
                      indices$within_id2_char[i]]
    indices$obseff[i] <- o[indices$within_id1_char[i],
                           indices$within_id2_char[i]]
  }

  if (sum(indices$b) != sum(unlist(lapply(bmats, sum)))) {
    stop("error at processing interaction data from input matrices")
  }

  # indi table
  inditable <- indilist
  inditable$trueindisoc <- NA

  if (forsimeval) {
    if (is.null(dyaddata)) stop("require list of true dyad data to assemble evaluation data")
    if (is.null(indidata)) stop("require list of true indi data to assemble evaluation data")
    indices$truedyadsoc <- NA

    for (i in seq_along(groups)) {
      g <- names(groups)[i]
      sel <- which(indices$period_char == g)
      indices$truedyadsoc[sel] <- dyaddata[[g]]$dyad_soc_vals
      sel <- which(inditable$period_char == g)
      inditable$trueindisoc[sel] <- indidata[[g]]$indi_soc_vals
    }
    return(list(indices = indices, inditable = inditable))
  }

  d <- list(n_pwise_dyads = nrow(indices),
            n_pwise_ids = nrow(inditable),
            n_periods = length(groups),
            behav_types = btyp_num,
            n_beh = 1,
            interactions = matrix(indices$b),
            n_ids_perperiod = as.numeric(table(inditable$period)),
            n_dyads_perperiod = as.numeric(table(indices$period)),
            obseff = matrix(indices$obseff),
            obseff_int = matrix(as.integer(indices$obseff)),
            index_period = indices$period,
            index_period_individual = inditable$period,
            index_dyad = indices$dyad,
            index_id1 = indices$id1,
            index_id2 = indices$id2,
            prior_matrix = matrix(ncol = 2, nrow = length(groups)),
            n_cors = 0,
            is_multi_manygroups = 1
  )
  names(d$index_period) <- indices$period_char
  names(d$index_id1) <- indices$id1_char
  names(d$index_id2) <- indices$id2_char
  names(d$index_dyad) <- indices$dyad_char
  names(d$index_period_individual) <- inditable$id_char
  names(d$n_dyads_perperiod) <- names(bmats)
  names(d$n_ids_perperiod) <- names(bmats)

  for (i in seq_along(groups)) {
    d$prior_matrix[i, ] <- make_prior(
      response = indices$b[indices$period == i],
      type = "count",
      obseff = indices$obseff[indices$period == i]
    )
  }
  # lapply(d, head)
  d
}
