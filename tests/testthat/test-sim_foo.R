test_that("function works", {
  res <- sim_foo(n_ids = 6, silent = FALSE, n_beh = 1)
  expect_vector(res)
  res <- sim_foo(n_ids = 6, silent = TRUE, n_beh = 1)
  expect_vector(res)
})
