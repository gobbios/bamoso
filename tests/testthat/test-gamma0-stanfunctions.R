f <- system.file("extdata/interaction_model.stan", package = "bamoso")
td <- tempdir()
mod <- cmdstan_model(stan_file = f, compile = TRUE, dir = td)
mod$expose_functions(global = TRUE)







