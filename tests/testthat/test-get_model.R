library(bamoso)

test_that("get_model return a CmdStanModel", {
  m <- get_model(type = "simple")
  expect_s3_class(m, "CmdStanModel")
  m <- get_model(type = "indi_cat")
  expect_s3_class(m, "CmdStanModel")
})
