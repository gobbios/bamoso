# stan data is generated currently in three functions:
# - generate_data
# - make_stan_data_from_associations
# - make_stan_data_from_matrices

# testing here whether the necessary elements are produced in each function


# no confounders -----------------------
behav_types <- sample(c("count", "prop", "dur_gamma", "dur_beta"), 2)
x <- generate_data(n_beh = 2, n_ids = 6, beh_intercepts = runif(2, -1, 1),
                   behav_types = behav_types)
mats <- x$processed$interaction_matrices

s1 <- x$standat
s2 <- make_stan_data_from_matrices(mats = mats, behav_types = behav_types)

res <- data.frame(list_item = names(s1), is_there = FALSE, is_equal = NA)

for (i in 1:nrow(res)) {
  if (res$list_item[i] %in% names(s2)) {
    res$is_there[i] <- TRUE
    res$is_equal[i] <- identical(s1[[res$list_item[i]]], s2[[res$list_item[i]]])
  }
}


test_that("data items are present in both simulated output and make_stan_data_from_matrices", {
  expect_true(all(res$is_there))
})

x <- sim_data_whitehead2015(n_ind = 5, n_periods = 9, ignore_temporal = TRUE)
assos <- x$assolist
s3 <- make_stan_data_from_association(assos)

res <- data.frame(list_item = names(s1), is_there = FALSE)
for (i in 1:nrow(res)) {
  if (res$list_item[i] %in% names(s3)) {
    res$is_there[i] <- TRUE
  }
}

test_that("data items are present in both simulated output and make_stan_data_from_association", {
  expect_true(all(res$is_there))
})



# with individual-level confounder -----------------------
behav_types <- sample(c("count", "prop", "dur_gamma", "dur_beta"), 2)
x <- generate_data(n_beh = 2, n_ids = 6, beh_intercepts = runif(2, -1, 1),
                   behav_types = behav_types, indi_cat_slope = 1)
mats <- x$processed$interaction_matrices

s1 <- x$standat
s2 <- make_stan_data_from_matrices(mats = mats, behav_types = behav_types, indi_cat_pred = x$input_data$indi_data$feature_cat)

res <- data.frame(list_item = names(s1), is_there = FALSE, is_equal = NA)

for (i in 1:nrow(res)) {
  if (res$list_item[i] %in% names(s2)) {
    res$is_there[i] <- TRUE
    res$is_equal[i] <- identical(s1[[res$list_item[i]]], s2[[res$list_item[i]]])
  }
}


test_that("data items are present in both simulated output and make_stan_data_from_matrices", {
  expect_true(all(res$is_there))
})

x <- sim_data_whitehead2015(n_ind = 5, n_periods = 9, ignore_temporal = TRUE)
assos <- x$assolist
s3 <- make_stan_data_from_association(assos, indi_cat_pred = sample(c(0, 1), ncol(assos), TRUE))

res <- data.frame(list_item = names(s1), is_there = FALSE)
for (i in 1:nrow(res)) {
  if (res$list_item[i] %in% names(s3)) {
    res$is_there[i] <- TRUE
  }
}

test_that("data items are present in both simulated output and make_stan_data_from_association", {
  expect_true(all(res$is_there))
})
