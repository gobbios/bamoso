#' generate simulated affiliation data
#'
#' @param n_ids numeric, number of individuals
#' @param n_beh numeric, number of behaviors
#' @param behav_types character vector of length \code{n_beh} that describes
#'        the kind of data to be generated. Possible values are \code{"prop"},
#'        \code{"count"}, \code{"dur_gamma"} and \code{"dur_beta"}.
#' @param indi_sd numeric, the SD for the individual component. Must be
#'                positive. Default is a random value
#'                (\code{runif(1, 0, 4)}).
#' @param dyad_sd numeric, the SD for the dyadic component. Must be
#'                positive. Default is a random value
#'                (\code{runif(1, 0, 4)}).
#' @param beh_intercepts numeric of the same length as \code{n_beh}: the
#'                       intercepts on the linear scale for each behavior
#'                       (the group-level average). Default is random
#'                       (\code{rnorm(n_beh, -1, 2)}).
#' @param dyadic_covariate_slope numeric, default is 0. Fits an additional
#'        covariate (stdnormal) with this slope to contribute to the linear
#'        predictor.
#' @param prop_trials numeric of length 1 or 2, default is \code{100}, i.e.
#'        per dyad there are 100 'trials' observed. If numeric of length 2,
#'        each dyad gets its own number of trials ranging between the two values
#'        supplied. Only relevant for data generated as proportions.
#' @param disp_pars_gamma numeric of the same length as \code{n_beh}.
#'                        Dispersion parameter(s) for any behavior with
#'                        \code{behav_types = "dur_gamma"}. See details.
#' @param disp_pars_beta numeric of the same length as \code{n_beh}.
#'                       Dispersion parameter(s) for any behavior with
#'                       \code{behav_types = "dur_beta"}. See details.
#' @param count_obseff numeric of length 1 or 2, default is \code{1}, i.e.
#'        per dyad there is 1 unit of observation effort (think 'hours').
#'        If numeric of length 2, each dyad get its own observation effort,
#'        sampled as uniform real between the two values. This is only relevant
#'        for data generated as counts.
#' @param exact logical, default is \code{TRUE}: should the varying intercepts
#'        for \code{indi_sd} and \code{dyad_sd} be rescaled so that they
#'        have means of 0 and exact SD as supplied.
#'
#' @details Currently, four data types/distributions via \code{behav_types}
#'          are supported are \code{"count"}, \code{"prop"},
#'          \code{"dur_gamma"} and \code{"dur_beta"}.
#'
#'          For the dispersion parameters, the input vector must be of the same
#'          length and the indexing must match \code{behav_types}. For example,
#'          \code{behav_types = c("count", "dur_gamma")} requires a vector of
#'          length 2, e.g. \code{disp_pars_gamma = c(0, 0.6)}, where the
#'          first item will be ignored, i.e. only the second entry is
#'          relevant.
#'
#' @importFrom stats cor density rnorm rpois runif rbinom rgamma
#' @importFrom stats plogis rbeta
#' @importFrom utils combn
#' @importFrom Rdpack reprompt
#'
#' @return a list with three named items, which are lists too:
#'  \itemize{
#'    \item{\code{$processed}:} {misc output generated, most importantly
#'                              interactions in matrix format}
#'    \item{\code{$standat}:} {standata to be passed to
#'                            \code{\link{sociality_model}} (see
#'                            \code{\link{make_stan_data_from_matrices}})}
#'    \item{\code{$input_data}:} {the input data supplied via arguments or
#'                               generated randomly in the function if using
#'                               the arguments' defaults}
#'  }
#'
#' @export
#' @examples
#' ex <- generate_data(n_ids = 7, n_beh = 2,
#'                     behav_types = c("count", "dur_gamma"),
#'                     indi_sd = 1.2, dyad_sd = 0.8,
#'                     disp_pars_gamma = c(0, 0.6),
#'                     beh_intercepts = c(1.4, -0.7))
#' # two interaction matrices:
#' ex$processed$interaction_matrices[[1]]
#' ex$processed$interaction_matrices[[2]]
#'

# n_ids=8
# indi_sd=2.5
# dyad_sd=1.2
# n_beh=2
# beh_intercepts=c(1.2, -2)
# exact=TRUE

generate_data <- function(n_ids = NULL,
                          n_beh = NULL,
                          behav_types = "count",
                          indi_sd = NULL,
                          dyad_sd = NULL,
                          dyadic_covariate_slope = 0, # slope for additional dyadic nuissance parameter
                          disp_pars_gamma = NULL, # for durations as gamma
                          disp_pars_beta = NULL, # for durations as beta
                          beh_intercepts = NULL,
                          prop_trials = 100, # number of trials for binomial (proportion data)
                          count_obseff = 1, # observation effort for count data
                          exact = TRUE) {

  # those used to be function arguments, but were removed
  indi_intercept <- 0
  dyad_intercept <- 0

  if (is.null(n_ids)) n_ids <- sample(5:30, 1)
  if (is.null(n_beh)) n_beh <- sample(1:4, 1)
  if (is.null(indi_sd)) indi_sd <- runif(1, 0, 4)
  if (is.null(dyad_sd)) dyad_sd <- runif(1, 0, 4)
  if (is.null(disp_pars_gamma)) disp_pars_gamma <- runif(n_beh, 0, 1)
  if (is.null(disp_pars_beta)) disp_pars_beta <- runif(n_beh, 5, 30)
  if (is.null(beh_intercepts)) beh_intercepts <- rnorm(n_beh, -1, 2)

  # prelims
  n_dyads <- n_ids * (n_ids - 1) / 2
  dyads <- t(combn(seq_len(n_ids), 2))

  # individual sociality
  indi_soc_vals <- rnorm(n = n_ids, mean = indi_intercept, sd = indi_sd)
  # dyadic sociality
  dyad_soc_vals <- rnorm(n = n_dyads, mean = dyad_intercept, sd = dyad_sd)


  # select behavior scale types, currently:
  #  - count (as in poisson)
  #  - prop (as in proportion): still technically a count (out of 100 'trials')
  #  - duration (not yet implemented...)

  if (length(behav_types) == 1) {
    if (behav_types == "random") {
      btypes <- sample(c("count", "prop", "dur_gamma"), n_beh, replace = TRUE)
    } else {
      btypes <- rep(behav_types, n_beh)
    }
  }

  if (length(behav_types) > 1) {
    if (length(behav_types) != n_beh) stop("number of behavior types does not match number of behaviors")
    if (any(behav_types == "random")) stop("can't have random behavior type with more than one behavior (need to specify for each behavior individually)")
    btypes <- behav_types
  }

  # and prep for output (convert to numeric)
  behav_types_num <- rep(1, n_beh) # counts (pois as default)
  behav_types_num[btypes == "prop"] <- 2
  behav_types_num[btypes == "dur_gamma"] <- 3
  behav_types_num[btypes == "dur_beta"] <- 4

  # indexing for optional shape/dispersion parameters
  gamma_shape <- logical(n_beh)
  gamma_shape[btypes == "dur_gamma"] <- TRUE
  gamma_shape_pos <- numeric(n_beh)

  beta_shape <- logical(n_beh)
  beta_shape[btypes == "dur_beta"] <- TRUE
  beta_shape_pos <- numeric(n_beh)




  # temporary observation effort
  obseff <- matrix(ncol = n_beh, nrow = n_dyads, 1)
  for (i in seq_len(n_beh)) {
    if (btypes[i] == "prop") {
      if (length(prop_trials) == 1) {
        obseff[, i] <- prop_trials
      } else {
        obseff[, i] <- sample(as.integer(prop_trials[1]):as.integer(prop_trials[2]), size = n_dyads, replace = TRUE)
      }
    }
    if (btypes[i] == "count") {
      if (length(count_obseff) == 1) {
        obseff[, i] <- count_obseff
      } else {
        obseff[, i] <- runif(n = n_dyads,
                             min = count_obseff[1],
                             max = count_obseff[2])
      }
    }
  }


  # random dyadic covariate (think: relatedness)
  dyadic_covariate_predictor <- as.numeric(scale(rnorm(n = n_dyads,
                                                       mean = 0,
                                                       sd = 1)))

  # exact distributions
  if (exact) {
    indi_soc_vals <- as.numeric(scale(indi_soc_vals) * indi_sd + indi_intercept)
    dyad_soc_vals <- as.numeric(scale(dyad_soc_vals) * dyad_sd + dyad_intercept)
  }

  lp <- indi_intercept + dyad_intercept +
    0.5 * (indi_soc_vals[dyads[, 1]] + indi_soc_vals[dyads[, 2]]) +
    dyad_soc_vals +
    dyadic_covariate_slope * dyadic_covariate_predictor

  interactions <- matrix(ncol = n_beh, nrow = n_dyads)
  for (i in seq_len(n_beh)) {
    lp_b <- lp + beh_intercepts[i]
    if (btypes[i] == "count") {
      interactions[, i] <- rpois(n = n_dyads,
                                 lambda = exp(lp_b + log(obseff[, i])))
    }
    if (btypes[i] == "prop") {
      interactions[, i] <- rbinom(n = n_dyads,
                                  size = obseff[, i],
                                  prob = exp(lp_b) / (1 + exp(lp_b))) # inverse logit...
    }
    if (btypes[i] == "dur_gamma") {
      # first translate mean (and variance) into shape and rate
      dispersion <- 1 / disp_pars_gamma[i]
      lp_offs <- exp(lp_b + log(obseff[, i]))
      variance <- dispersion * (lp_offs^2)
      shape <- (lp_offs^2) / variance
      rate <- lp_offs / variance
      interactions[, i] <- rgamma(n = n_dyads, shape = shape, rate = rate)

      gamma_shape_pos[i] <- i
    }

    if (btypes[i] == "dur_beta") {
      # first translate mean (and variance) into shape1 and shape2 (for rbeta())
      lp_offs <- lp_b
      # mean needs to be on logit scale
      mu <- plogis(lp_offs)
      phi <- disp_pars_beta[i]

      interactions[, i] <- rbeta(n = n_dyads,
                                 shape1 = mu * phi,
                                 shape2 = (1 - mu) * phi)
      # make sure no 0 and 1 in the results
      # interactions[, i] <- (interactions[, i] * (n_dyads - 1) + 0.5) / n_dyads

      beta_shape_pos[i] <- i
    }


  }

  # create default priors for behaviors
  # at this point 'interactions' is not yet split into discrete vs continuous
  prior_matrix <- matrix(ncol = 2, nrow = n_beh)
  for (i in 1:n_beh) {
    response <- interactions[, i]
    prior_matrix[i, ] <- make_prior(response = response,
                                    type = behav_types[i],
                                    obseff = obseff[, i])
  }

  # dsi
  aux <- interactions / obseff
  aux <- apply(aux, 2, function(y) y / mean(y))
  dsi <- rowMeans(aux)

  # sum of top 3
  sum_top3 <- numeric(n_ids)
  for (i in seq_len(n_ids)) {
    x <- dsi[apply(dyads, 1, function(x) any(x == i))]
    sum_top3[i] <- sum(sort(x, decreasing = TRUE)[1:3])
  }

  # interaction matrices
  imats <- vector("list", n_beh)
  for (i in seq_len(n_beh)) {
    m <- matrix(ncol = n_ids, nrow = n_ids, 0)
    for (k in seq_len(nrow(dyads))) {
      m[dyads[k, 1], dyads[k, 2]] <- interactions[k, i]
    }
    imats[[i]] <- m
  }

  # flag for when there is an empty interaction matrix generated
  empty <- any(unlist(lapply(imats, sum)) == 0)
  if (empty) {
    if (interactive()) {
      warning("there was an empty matrix generated",
              call. = FALSE)
    } else {
      message("there was an empty matrix generated")
    }
  }
  #
  gamma_shape_pos_mod <- numeric(n_beh)
  gamma_shape_pos_mod[gamma_shape_pos > 0] <- seq_len(sum(gamma_shape_pos > 0))
  beta_shape_pos_mod <- numeric(n_beh)
  beta_shape_pos_mod[beta_shape_pos > 0] <- seq_len(sum(beta_shape_pos > 0))

  # two tweaks:
  #   behav_types needs to be of at least length 2 (so it can be declared
  #               a non-scalar for stan): hence a 0 is appended
  #   obseff_int needs to be declared as an integer matrix so that it
  #              can account for binomial trials

  # create placeholder id names
  id_codes <- as.character(seq_len(n_ids))
  vec <- rep(0, n_ids)
  names(vec) <- id_codes

  # code behave_types as named vector
  bdata <- c(behav_types_num, 0)
  names(bdata) <- c(behav_types, "0")

  standat <- list(id1 = dyads[, 1],
                  id2 = dyads[, 2],
                  behav_types = bdata,
                  interactions = apply(interactions, 2, as.integer),
                  interactions_cont = interactions,
                  n_beh = n_beh,
                  n_ids = n_ids,
                  n_dyads = n_dyads,
                  dyads_navi = as.matrix(dyads[, 1:2]),
                  obseff = obseff,
                  prior_matrix = prior_matrix,
                  obseff_int = apply(obseff, 2, as.integer),
                  gamma_shape_pos = gamma_shape_pos_mod,
                  gamma_shape_n = sum(gamma_shape_pos > 0),
                  beta_shape_pos = beta_shape_pos_mod,
                  beta_shape_n = sum(beta_shape_pos > 0),
                  id_codes = vec,
                  beh_names = rep(0, n_beh)
                  )

  names(standat$beh_names) <- paste0("behav_", LETTERS[seq_len(n_beh)])

  list(standat = standat,
       input_data = list(behav_type = btypes,
                         empty = empty,
                         indi_soc_vals = indi_soc_vals,
                         indi_sums = 0.5 * (indi_soc_vals[dyads[, 1]] +
                                              indi_soc_vals[dyads[, 2]]),
                         dyad_soc_vals = dyad_soc_vals,
                         dyads = dyads,
                         indi_sd = indi_sd,
                         indi_intercept = indi_intercept,
                         dyad_sd = dyad_sd,
                         dyad_intercept = dyad_intercept,
                         dyadic_covariate_predictor = dyadic_covariate_predictor,
                         dyadic_covariate_slope = dyadic_covariate_slope,
                         disp_pars_gamma = disp_pars_gamma,
                         disp_pars_beta = disp_pars_beta,
                         beh_intercepts = beh_intercepts,
                         obseff = obseff),
       processed = list(interaction_cor = ifelse(all(colSums(interactions) > 0),
                                                 cor(interactions, method = "s"),
                                                 NA),
                        dsi = dsi,
                        sum_top3 = sum_top3,
                        interaction_matrices = imats))

}
