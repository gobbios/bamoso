test_that("multiple datasets/groups/periods", {
  n <- 3
  groups <- bamoso:::generate_groups(n_groups = n, min_nids = 20, max_nids = 25)
  (indisds <- round(exp(rnorm(n, -1, 0.7)), 1))
  (dyadsds <- round(exp(rnorm(n, -0.5, 0.6)), 1))
  (bintercepts <- round(rnorm(n, -2, 0.8), 1))

  xx <- bamoso:::generate_data_multi(groups, indisds, dyadsds, bintercepts)

  d <- bamoso:::make_stan_data_from_matrices_multi(
    bmats = xx$bmats,
    omats = xx$omats,
    btypes = "count"
  )

  r <- sociality_model(d, parallel_chains = 4, refresh = 500, show_exceptions = FALSE)

  expect_error(pp_model(r, group = "ABC"))
  expect_no_error(pp_model(r, group = "P01"))
  expect_no_error(pp_model(r, group = c("P03", "P01")))
  expect_no_error(pp_model(r, group = c("P02")))
  expect_no_error(pp_model_stat(r, stat = "max"))
  expect_no_error(pp_model_stat(r, stat = "max", group = "P01"))
  expect_no_error(pp_model_stat(r, stat = "max", group = c("P03", "P01")))
  expect_error(pp_model_stat(r, stat = "max", group = "ABC"))

})

test_that("multiple datasets/groups/periods", {
  n <- sample(2:4, 1)
  groups <- bamoso:::generate_groups(n_groups = n, min_nids = 5, max_nids = 7)
  (indisds <- round(exp(rnorm(n, -1, 0.7)), 1))
  (dyadsds <- round(exp(rnorm(n, -0.5, 0.6)), 1))
  (bintercepts <- round(rnorm(n, -2, 0.8), 1))

  xx <- bamoso:::generate_data_multi(groups, indisds, dyadsds, bintercepts)

  d <- bamoso:::make_stan_data_from_matrices_multi(
    bmats = xx$bmats,
    omats = xx$omats,
    btypes = "count"
  )

  p0 <- list(indi_sd = 2, dyad_sd = 2)
  expect_no_error(sociality_model(d, priors = p0, parallel_chains = 2, chains = 2, iter_warmup = 300, iter_sampling = 100, refresh = 0, show_exceptions = FALSE))
  p1 <- list(indi_sd = runif(n + 1, 0, 1))
  expect_error(sociality_model(d, priors = p1, parallel_chains = 2, chains = 2, iter_warmup = 300, iter_sampling = 100, refresh = 0, show_exceptions = FALSE))
  p2 <- list(indi_sd = runif(n, 0, 1))
  expect_error(sociality_model(d, priors = p2, parallel_chains = 2, chains = 2, iter_warmup = 300, iter_sampling = 100, refresh = 0, show_exceptions = FALSE))
  p3 <- list(indi_sd = runif(n, 1, 2), dyad_sd = runif(n, 1, 2))
  expect_no_error(sociality_model(d, priors = p3, parallel_chains = 2, chains = 2, iter_warmup = 300, iter_sampling = 100, refresh = 0, show_exceptions = FALSE))

})

