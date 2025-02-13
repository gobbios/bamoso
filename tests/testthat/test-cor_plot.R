test_that("cor_plot doesn't crash", {
  capture.output(suppressMessages(model <- fit_example_model(1)))
  expect_null(cor_plot(model, asp = 1, xlab = "A", ylab = "B"))
  expect_null(cor_plot(model, greg = FALSE, xlim = c(-3, 3), asp = 1))
})
