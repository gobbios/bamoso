#' wrapper to provide model executable
#'
#' @return a cmdstan model object
#' @export
#'
#' @importFrom cmdstanr cmdstan_model
#' @importFrom Rdpack reprompt

get_model <- function() {
  f <- system.file("extdata/interaction_model.stan", package = "bamoso")
  mod <- cmdstan_model(stan_file = f, compile = TRUE, dir = NULL)
  mod
}
