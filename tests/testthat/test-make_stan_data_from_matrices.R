
# occasionally, an empty matrix is produced (mostly in the gamma setting)
good_to_go <- FALSE
while (!good_to_go) {
  x <- generate_data(n_beh = 2, dyad_sd = 1, indi_sd = 1, beh_intercepts = runif(2, -1, 1),
                     behav_types = sample(c("count", "prop", "dur_gamma", "dur_beta"), 2))
  if (!x$input_data$empty) good_to_go <- TRUE
}

x$input_data$empty
x$input_data$behav_type
bmats <- list(b1 = x$processed$interaction_matrices[[1]],
              b2 = x$processed$interaction_matrices[[2]])

standat <- make_stan_data_from_matrices(mats = bmats,
                                        behav_types = x$input_data$behav_type)

test_that("make_stan_data_from_matrices does not produce NA values", {
  expect_true(all(unlist(lapply(standat, function(x) !any(is.na(x))))))
})

test_that("make_stan_data_from_matrices contains the correct list elements", {
  expect_true("id1" %in% names(standat))
  expect_true("id2" %in% names(standat))
  expect_true("interactions" %in% names(standat))
  expect_true("interactions_cont" %in% names(standat))
  expect_true("n_beh" %in% names(standat))
  expect_true("n_ids" %in% names(standat))
  expect_true("n_dyads" %in% names(standat))
  expect_true("dyads_navi" %in% names(standat))
  expect_true("gamma_shape_pos" %in% names(standat))
  expect_true("gamma_shape_n" %in% names(standat))
  expect_true("beta_shape_pos" %in% names(standat))
  expect_true("beta_shape_n" %in% names(standat))
  expect_true("obseff" %in% names(standat))
  expect_true("obseff_int" %in% names(standat))
  expect_true("prior_matrix" %in% names(standat))
  expect_true("id_codes" %in% names(standat))
  expect_true("beh_names" %in% names(standat))
  expect_true("behav_types" %in% names(standat))
})


# dimensions match
nrow_interactions <- nrow(standat$interactions)
nrow_interactions_cont <- nrow(standat$interactions_cont)

