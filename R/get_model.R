#' wrapper to provide model executable
#'
#' @param type character, specify which model
#' @return a cmdstan model object
#' @export
#'
#' @importFrom cmdstanr cmdstan_model
#' @importFrom Rdpack reprompt

get_model <- function(type = c("simple", "indi_cat")) {
  if (type == "simple") {
    f <- system.file("extdata/interaction_model.stan", package = "bamoso")
    mod <- cmdstan_model(stan_file = f, compile = TRUE, dir = NULL)
  }
  if (type == "indi_cat") {
    f <- system.file("extdata/interaction_model_indicat.stan", package = "bamoso")
    mod <- cmdstan_model(stan_file = f, compile = TRUE, dir = NULL)
  }

  mod
}
