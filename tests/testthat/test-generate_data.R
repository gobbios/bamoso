# try to consistently produce an empty matrix which should trigger a warning (if interactive)

test_that("warning/message is generated if any matrix is empty", {
  if (interactive()) {
    expect_warning(generate_data(n_beh = 4, n_ids = 6, dyad_sd = 1, indi_sd = 1, beh_intercepts = runif(4, -100, -90),
                                 behav_types = c("count", "prop", "dur_gamma", "dur_beta")))
  } else {
    expect_message(generate_data(n_beh = 4, n_ids = 6, dyad_sd = 1, indi_sd = 1, beh_intercepts = runif(4, -100, -90),
                                 behav_types = c("count", "prop", "dur_gamma", "dur_beta")))
  }

})
