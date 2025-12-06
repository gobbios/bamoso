

x <- generate_data(n_ids = 12, n_beh = 2, behav_types = c("prop", "count"),
                   indi_sd = 1, dyad_sd = 1, indi_covariate_slope = 0.2,
                   indi_cat_slope = -0.5, dyad_covariate_slope = 0,
                   dyad_cat_slope = 0, beh_intercepts = c(0, 0),
                   prop_trials = 500, count_obseff = 10)
b <- list(b1 = x$processed$interaction_matrices[[1]],
          b2 = x$processed$interaction_matrices[[2]])
o <- x$processed$obseff_matrices

# basic model
d1 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"),
                                   obseff = o)
f1 <- sociality_model(d1, parallel_chains = 4, iter_warmup = 500,
                      iter_sampling = 200, show_exceptions = FALSE,
                      show_messages = FALSE, diagnostics = NULL)
expect_null(suppressMessages(summary(f1)))

# corr model
d2 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"),
                                   obseff = o, correlations = TRUE)
f2 <- sociality_model(d2, parallel_chains = 4, iter_warmup = 500,
                      iter_sampling = 200, show_exceptions = FALSE,
                      show_messages = FALSE, diagnostics = NULL)
expect_null(suppressMessages(summary(f2)))

# sans dyadic model
d3 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"),
                                   obseff = o)
f3 <- sociality_model(d3, parallel_chains = 4, iter_warmup = 500,
                      iter_sampling = 200, show_exceptions = FALSE,
                      show_messages = FALSE, diagnostics = NULL,
                      sans_dyadic = TRUE)
expect_null(suppressMessages(summary(f3)))

# covar model
d4 <- make_stan_data_from_matrices(
  mats = b,
  behav_types = c("prop", "count"),
  obseff = o,
  indi_cat_pred = x$input_data$indi_data$feature_cat
  )
f4 <- sociality_model(d4, parallel_chains = 4, iter_warmup = 500,
                      iter_sampling = 200, show_exceptions = FALSE,
                      show_messages = FALSE, diagnostics = NULL)

# grooming model over four groups
gmod <- fit_example_model(num = "grooming1")

test_that("model type is recognized after fitting", {
  expect_true(f1$modeltype == "simple")
  expect_true(f2$modeltype == "cor_mod")
  expect_true(f3$modeltype == "sans_dyadic")
  expect_true(f4$modeltype == "with_preds")
  expect_true(gmod$modeltype == "multi_manygroups")
})


test_that("summary doesn't throw errors", {
  zz <- capture.output(expect_null(suppressMessages(summary(f1))))
  zz <- capture.output(expect_null(suppressMessages(summary(f2))))
  zz <- capture.output(expect_null(suppressMessages(summary(f3))))
  zz <- capture.output(expect_null(suppressMessages(summary(f4))))
  zz <- capture.output(expect_no_error(suppressWarnings(summary(gmod))))
})
