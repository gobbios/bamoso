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
