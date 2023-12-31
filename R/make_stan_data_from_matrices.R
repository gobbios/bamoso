#' convert interaction matrices into data list ready for stan
#'
#' @param mats list with (named) square interaction matrices
#' @param obseff numeric vector with observation effort. Must be positive
#'               and one value for each individual, e.g. focal observation
#'               time or number of samples. Can also be a square matrix in which
#'               dyadic observation time is already coded for each dyad.
#'               Can also be a list of such vectors or matrices, which then
#'               is mapped to each behavior (in that case this list needs to be
#'               of the same length as \code{mats}).
#' @param behav_types character vector of length \code{n_beh} that describes
#'                    the kind of data. Possible values are \code{"prop"},
#'                    \code{"count"}, \code{"dur_gamma"} and \code{"dur_beta"}.
#'                    At its default all behaviors are considered
#'                    to be \code{"count"}.
#'
#' @details If supplied via column names in \code{mats[[1]]}, id codes are
#'          present in the output as names of vector of zeros.
#'
#'
#' @importFrom cmdstanr cmdstan_model
#' @return a named list that can be handed over to Stan as \code{data} via
#'         \code{\link{sociality_model}}:
#'  \itemize{
#'     \item{\code{$id1}} {integer index of individual 1}
#'     \item{\code{$id2}} {integer index of individual 2}
#'     \item{\code{$interactions}} {matrix of interactions (integer)}
#'     \item{\code{$interactions_cont}} {matrix of interactions (continuous)}
#'     \item{\code{$n_beh}} {integer number of behaviors}
#'     \item{\code{$n_ids}} {integer number of individuals}
#'     \item{\code{$n_dyads}} {integer number of dyads}
#'     \item{\code{$dyads_navi}} {integer index dyads}
#'     \item{\code{$obseff}} {matrix with continuous observation effort}
#'     \item{\code{$obseff_int}} {matrix with integer observation effort}
#'     \item{\code{$gamma_shape_n}} {number of response vectors with gamma
#'                                       likelihood}
#'     \item{\code{$gamma_shape_pos}} {integer index for responses with gamma
#'                                         likelihood}
#'     \item{\code{$beta_shape_n}} {number of response vectors with beta
#'                                      likelihood}
#'     \item{\code{$beta_shape_pos}} {integer index for responses with
#'                                        beta likelihood}
#'     \item{\code{$prior_matrix}} {matrix with priors for intercepts}
#'     \item{\code{$generate_predictions}} {currently unused}
#'     \item{\code{$id_codes}} {vector of 0's with names corresponding to
#'                                  individuals' codes}
#'     \item{\code{$beh_names}} {vector of 0's with names corresponding to
#'                                   behavior codes}
#'  }
#' @export
#'
#' @examples
#' mats <- vector(mode = "list", length = 2)
#' m1 <- matrix(rpois(64, 1), ncol = 8)
#' diag(m1) <- 0
#' m2 <- matrix(rpois(64, 5), ncol = 8)
#' diag(m2) <- 0
#' colnames(m1) <- letters[1:8]
#'
#' matlist <- list(beh1 = m1, beh2 = m2)
#' make_stan_data_from_matrices(mats = matlist)
#'
#' o <- rep(0.5, 8)
#' o[1:2] <- c(1.5, 10.5)
#' make_stan_data_from_matrices(mats = matlist, obseff = o)$obseff

make_stan_data_from_matrices <- function(mats,
                                         behav_types = NULL,
                                         obseff = NULL) {
  n_ids <- ncol(mats[[1]])
  n_dyads <- n_ids * (n_ids - 1) * 0.5
  n_beh <- length(mats)

  # deal with behavior types
  if (is.null(behav_types)) {
    behav_types <- rep("count", n_beh)
  }
  if (length(behav_types) != n_beh) {
    stop("require exactly ", n_beh, " 'behave_types' specified")
  }
  if (!all(behav_types %in% c("count", "prop", "dur_gamma", "dur_beta"))) {
    stop("unknown behavior type supplied")
  }

  # and prep for output (convert to numeric)
  behav_types_num <- rep(1, n_beh) # counts (pois as default)
  behav_types_num[behav_types == "prop"] <- 2
  behav_types_num[behav_types == "dur_gamma"] <- 3
  behav_types_num[behav_types == "dur_beta"] <- 4

  # indexing for optional shape/dispersion parameters
  gamma_shape <- logical(n_beh)
  gamma_shape[behav_types == "dur_gamma"] <- TRUE
  gamma_shape_pos <- numeric(n_beh)

  beta_shape <- logical(n_beh)
  beta_shape[behav_types == "dur_beta"] <- TRUE
  beta_shape_pos <- numeric(n_beh)

  for (i in seq_len(n_beh)) {
    if (behav_types[i] == "dur_gamma") {
      gamma_shape_pos[i] <- i
    }
    if (behav_types[i] == "dur_beta") {
      beta_shape_pos[i] <- i
    }
  }

  gamma_shape_pos_mod <- numeric(n_beh)
  gamma_shape_pos_mod[gamma_shape_pos > 0] <- seq_len(sum(gamma_shape_pos > 0))
  beta_shape_pos_mod <- numeric(n_beh)
  beta_shape_pos_mod[beta_shape_pos > 0] <- seq_len(sum(beta_shape_pos > 0))


  interactions <- matrix(ncol = length(mats), nrow = n_dyads)
  obseff_dat <- matrix(ncol = length(mats), nrow = n_dyads)
  index <- which(upper.tri(mats[[1]]), arr.ind = TRUE)
  if (is.null(obseff)) {
    sdata <- rep(0.5, n_ids)
  } else {
    sdata <- obseff
  }

  for (i in seq_along(mats)) {
    for (k in seq_len(n_dyads)) {
      interactions[k, i] <- mats[[i]][index[k, 1], index[k, 2]]
      interactions[k, i] <- interactions[k, i] + mats[[i]][index[k, 2],
                                                           index[k, 1]]
    }
  }

  # deal with observation effort -------------
  # 1) supplied as single vector or matrix (applies to all behaviors equally)
  if (!is.list(sdata)) {
    if (is.vector(sdata)) {
      # make sure order is the same in vector and interaction matrix columns
      if (!is.null(colnames(mats[[1]])) && !is.null(names(sdata))) {
        sdata <- sdata[colnames(mats[[1]])]
      }
      for (k in seq_len(n_dyads)) {
        obseff_dat[k, ] <- sdata[index[k, 1]] + sdata[index[k, 2]]
      }
    }
    if (is.matrix(sdata)) {
      obseff_dat[, ] <- sdata[index]
    }
  }

  # 2) supplied as list of vectors or matrices (mapped onto each behavior)
  if (is.list(sdata)) {
    if (length(sdata) != length(mats)) {
      stop("behaviour matrices and observation effort can't be matched")
    }
    for (i in seq_along(mats)) {
      if (is.vector(sdata[[i]])) {
        for (k in seq_len(n_dyads)) {
          obseff_dat[k, i] <- sdata[[i]][index[k, 1]] + sdata[[i]][index[k, 2]]
        }
      }
      if (is.matrix(sdata[[i]])) {
        obseff_dat[, i] <- sdata[[i]][index]
      }
    }
  }

  # set obseff for props to 100 if not supplied...
  if (is.null(obseff)) {
    for (i in seq_along(mats)) {
      if (behav_types[i] == "prop") {
        obseff_dat[, i] <- 100
      }
    }
  }

  # create default priors for behaviors
  # at this point 'interactions' is not yet split into discrete vs continuous
  prior_matrix <- matrix(ncol = 2, nrow = length(mats))
  for (i in seq_along(mats)) {
    response <- interactions[, i]
    prior_matrix[i, ] <- make_prior(response = response,
                                    type = behav_types[i],
                                    obseff = obseff_dat[, i])
  }

  if (!is.null(colnames(mats[[1]]))) {
    id_codes <- as.character(colnames(mats[[1]]))
  } else {
    id_codes <- as.character(seq_len(n_ids))
  }
  vec <- rep(0, n_ids)
  names(vec) <- id_codes

  if (!is.null(names(mats))) {
    beh_codes <- as.character(names(mats))
  } else {
    beh_codes <- paste0("behav_", LETTERS[seq_len(length(mats))])
  }
  vec_behaviors <- rep(0, length(beh_codes))
  names(vec_behaviors) <- beh_codes

  # code behave_types as named vector
  bdata <- c(behav_types_num, 0)
  names(bdata) <- c(behav_types, "0")

  list(id1 = index[, 1],
       id2 = index[, 2],
       behav_types = bdata,
       interactions = apply(interactions, 2, as.integer),
       interactions_cont = apply(interactions, 2, as.numeric),
       n_beh = length(mats),
       n_ids = n_ids,
       n_dyads = n_dyads,
       dyads_navi = as.matrix(index[, 1:2]),
       obseff = obseff_dat,
       obseff_int = apply(obseff_dat, 2, as.integer),
       gamma_shape_pos = gamma_shape_pos_mod,
       gamma_shape_n = sum(gamma_shape_pos > 0),
       beta_shape_pos = beta_shape_pos_mod,
       beta_shape_n = sum(beta_shape_pos > 0),
       prior_matrix = prior_matrix,
       id_codes = vec,
       beh_names = vec_behaviors)
}
