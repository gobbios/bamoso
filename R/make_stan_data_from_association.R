#' convert group-by-individual association data into data list ready for Stan
#'
#' @param asso_table a 0/1 matrix where columns correspond to individuals and
#'        rows to individual samples
#'
#' @details For now differences in sampling periods are ignored: each
#'          individual is assumed to be present (theoretically observable)
#'          in each sample.
#'
#'          If supplied via column names in \code{asso_table}, id codes are
#'          present in the output as names of vector of zeros.
#'
#' @return a list
#' @export
#'
#' @examples
#' # generate some data
#' x <- sim_data_whitehead2015(n_ind = 5, n_periods = 9, ignore_temporal = TRUE)
#' assos <- x$assolist
#' head(assos)
#'
#' make_stan_data_from_association(assos)
#'


make_stan_data_from_association <- function(asso_table) {
  # get dyadic counts
  ares <- asso_indices(asso_table)

  n_ids <- length(unique(c(ares[, "i1"], ares[, "i2"])))
  n_dyads <- n_ids * (n_ids - 1) * 0.5

  interactions <- matrix(ncol = 1, nrow = n_dyads)
  obseff_dat <- matrix(ncol = 1, nrow = n_dyads)
  m <- matrix(ncol = n_ids, nrow = n_ids)
  index <- which(upper.tri(m), arr.ind = TRUE)

  for (k in seq_len(n_dyads)) {
    xline <- index[k, 1] == ares[, "i1"] & index[k, 2] == ares[, "i2"]
    interactions[k, 1] <- ares[xline, "seen_together"]
    obseff_dat[k, 1] <- ares[xline, "seen_either"]
  }

  if (!is.null(colnames(asso_table))) {
    id_codes <- as.character(colnames(asso_table))
  } else {
    id_codes <- as.character(seq_len(n_ids))
  }
  if (all(id_codes == "")) {
    id_codes <- as.character(seq_len(n_ids))
  }
  vec <- rep(0, n_ids)
  names(vec) <- id_codes

  # make prior
  prior_matrix <- matrix(ncol = 2, nrow = 1)
  prior_matrix[1, ] <- make_prior(response = interactions[, 1],
                                  obseff = obseff_dat[, 1],
                                  type = "prop")
  vec_behaviors <- c(assoc = 0)

  # code behave_types as named vector
  bdata <- c(2, 0) # prop
  names(bdata) <- c("prop", "0")

  list(id1 = index[, 1],
       id2 = index[, 2],
       interactions = apply(interactions, 2, as.integer),
       interactions_cont = apply(interactions, 2, as.numeric),
       n_beh = 1,
       n_ids = n_ids,
       n_dyads = n_dyads,
       dyads_navi = as.matrix(index[, 1:2]),
       gamma_shape_pos = 0,
       gamma_shape_n = 0,
       beta_shape_pos = 0,
       beta_shape_n = 0,
       obseff = obseff_dat,
       obseff_int = apply(obseff_dat, 2, as.integer),
       prior_matrix = prior_matrix,
       id_codes = vec,
       beh_names = vec_behaviors,
       behav_types = bdata
       )
}
