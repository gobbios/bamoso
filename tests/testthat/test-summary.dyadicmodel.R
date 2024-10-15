

x <- generate_data(n_ids = 12, n_beh = 2, behav_types = c("prop", "count"), indi_sd = 1, dyad_sd = 1, indi_covariate_slope = 0.2, indi_cat_slope = -0.5, dyadic_covariate_slope = 0, dyadic_cat_slope = 0, beh_intercepts = c(0, 0), prop_trials = 500, count_obseff = 10)
b <- list(b1 = x$processed$interaction_matrices[[1]], b2 = x$processed$interaction_matrices[[2]])
o <- x$processed$obseff_matrices

# basic model
d1 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"), obseff = o)
f1 <- sociality_model(d1, parallel_chains = 4, iter_warmup = 500, iter_sampling = 200, show_exceptions = FALSE, show_messages = FALSE, diagnostics = NULL)
expect_null(suppressMessages(summary(f1)))

# corr model
d2 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"), obseff = o, correlations = TRUE)
f2 <- sociality_model(d2, parallel_chains = 4, iter_warmup = 500, iter_sampling = 200, show_exceptions = FALSE, show_messages = FALSE, diagnostics = NULL)
expect_null(suppressMessages(summary(f2)))

# sans dyadic model
d3 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"), obseff = o)
f3 <- sociality_model(d3, parallel_chains = 4, iter_warmup = 500, iter_sampling = 200, show_exceptions = FALSE, show_messages = FALSE, diagnostics = NULL, sans_dyadic = TRUE)
expect_null(suppressMessages(summary(f3)))

# covar model
d4 <- make_stan_data_from_matrices(mats = b, behav_types = c("prop", "count"), obseff = o, indi_cat_pred = x$input_data$indi_data$feature_cat)
f4 <- sociality_model(d4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 200, show_exceptions = FALSE, show_messages = FALSE, diagnostics = NULL)




test_that("summary doesn't throw errors", {
  expect_null(suppressMessages(summary(f1)))
  expect_null(suppressMessages(summary(f2)))
  expect_null(suppressMessages(summary(f3)))
  expect_null(suppressMessages(summary(f4)))
})
