library(bamoso)

test_that("get_model return a CmdStanModel", {
  m <- get_model()
  expect_s3_class(m, "CmdStanModel")
})
