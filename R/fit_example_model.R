
#' fit simple example models with simple example data
#'
#' @param num an integer
#' @importFrom utils data
#' @return a model object
#' @export
#'

fit_example_model <- function(num = 1) {
  if (num == "grooming1") {
    if (interactive()) {
      cat("macaque grooming model\n")
    }
    data("grooming", envir = environment()) # prevent loading in user workspace
    g <- list(ass = grooming$ass$groom,
              fus = grooming$fus$groom,
              nig = grooming$nig$groom,
              syl = grooming$syl$groom)
    o <- list(ass = grooming$ass$obseff,
              fus = grooming$fus$obseff,
              nig = grooming$nig$obseff,
              syl = grooming$syl$obseff)
    d <- bamoso:::make_stan_data_from_matrices_multi(
      bmats = g,
      btypes = "count",
      omats = o
      )
    r <- sociality_model(d,
                         parallel_chains = 4,
                         silent = TRUE,
                         seed = 3,
                         adapt_delta = 0.95,
                         iter_warmup = 400,
                         iter_sampling = 200,
                         refresh = 0,
                         show_exceptions = FALSE)
    return(r)

  }

  if (num == "grooming2") {
    if (interactive()) {
      cat("macaque grooming model with 6 fems\n")
    }
    data("grooming", envir = environment())
    g <- list(ass = grooming$ass$groom[1:6, 1:6],
              sylvanus = grooming$syl$groom[1:6, 1:6])
    o <- list(ass = grooming$ass$obseff[1:6, 1:6],
              sylvanus = grooming$syl$obseff[1:6, 1:6])
    d <- bamoso:::make_stan_data_from_matrices_multi(
      bmats = g,
      btypes = "count",
      omats = o
    )
    r <- sociality_model(d,
                         parallel_chains = 4,
                         silent = TRUE,
                         seed = 3,
                         adapt_delta = 0.95,
                         iter_warmup = 400,
                         iter_sampling = 200,
                         refresh = 0,
                         show_exceptions = FALSE)
    return(r)
  }

  if (num == 1) {
    if (interactive()) {
      cat("12 indis, 2 behaviors, cors generated, and fitted\n")
    }
    cors_indi <- matrix(c(0.7, 0.8, 0.8, 0.8), ncol = 2)
    cors_dyad <- matrix(c(0.7, 0.3, 0.3, 1.0), ncol = 2)

    x <- generate_data(n_ids = 12, n_beh = 2, beh_intercepts = c(0.5, 0.5),
                       behav_types = c("count", "count"),
                       indi_sd = cors_indi,
                       dyad_sd = cors_dyad)

    s <- x$standat

    r <- sociality_model(s, parallel_chains = 4, seed = 1, adapt_delta = 0.85,
                         show_exceptions = FALSE, refresh = 0,
                         iter_warmup = 500, iter_sampling = 100)
    return(r)
  }


  if (num == 2) {
    if (interactive()) {
      cat("12 indis, 2 behaviors, cors generated,",
          "but *not used* during data generation\n")
    }

    cors_indi <- matrix(c(0.7, 0.8, 0.8, 0.8), ncol = 2)
    cors_dyad <- matrix(c(0.7, 0.3, 0.3, 1.0), ncol = 2)

    x <- generate_data(n_ids = 12, n_beh = 2, beh_intercepts = c(0.5, 0.5),
                       behav_types = c("count", "count"),
                       indi_sd = cors_indi,
                       dyad_sd = cors_dyad)
    s <- x$standat
    s$n_cors <- 0

    r <- sociality_model(s, parallel_chains = 4, seed = 1, adapt_delta = 0.85,
                         show_exceptions = FALSE, refresh = 0,
                         iter_warmup = 500, iter_sampling = 100)
    return(r)
  }

  if (num == 3) {
    if (interactive()) {
      cat("12 indis, 3 behaviors, cors generated, and fitted\n")
    }

    cors_indi <- matrix(c(0.7, 0.8, 0.2,
                          0.8, 0.9, 0.4,
                          0.2, 0.4, 1.2),
                        ncol = 3)
    cors_dyad <- matrix(c(0.3, 0.8, -0.2,
                          0.8, 0.4, 0.4,
                          -0.2, 0.4, 1.0),
                        ncol = 3)

    x <- generate_data(n_ids = 12, n_beh = 3, beh_intercepts = c(0.5, 0.5, 0.5),
                       behav_types = c("count", "count", "count"),
                       indi_sd = cors_indi, dyad_sd = cors_dyad)

    s <- x$standat
    r <- sociality_model(s, parallel_chains = 4, seed = 1, adapt_delta = 0.85,
                         show_exceptions = FALSE, refresh = 0,
                         iter_warmup = 500, iter_sampling = 100)
    return(r)
  }

  if (num == 4) {
    if (interactive()) {
      cat("12 indis, 4 behaviors, cors generated, and fitted\n")
    }

    cors_indi <- matrix(c(1, 0.5, 0.5, -0.2,
                          0.5, 1, 0.5, 0.5,
                          0.5, 0.5, 1, 0.5,
                          -0.2, 0.5, 0.5, 1),
                        ncol = 4)
    cors_dyad <- matrix(c(1, 0.5, 0.5, -0.2,
                          0.5, 1, 0.5, 0.5,
                          0.5, 0.5, 1, 0.5,
                          -0.2, 0.5, 0.5, 1),
                        ncol = 4)

    x <- generate_data(n_ids = 12, n_beh = 4,
                       beh_intercepts = c(0.5, 0.5, 0.5, 0.5),
                       behav_types = c("count", "count", "count", "count"),
                       indi_sd = cors_indi, dyad_sd = cors_dyad)

    s <- x$standat
    r <- sociality_model(s, parallel_chains = 4, seed = 1, adapt_delta = 0.85,
                         show_exceptions = FALSE, refresh = 0,
                         iter_warmup = 500, iter_sampling = 100)
    return(r)
  }

  if (num == 5) {
    if (interactive()) {
      cat("12 indis, 2 behaviors, no cors\n")
    }

    x <- generate_data(n_ids = 12, n_beh = 2,
                       beh_intercepts = c(0.5, 0.5),
                       behav_types = c("count", "prop"),
                       indi_sd = 1.2, dyad_sd = 0.8)

    s <- x$standat
    r <- sociality_model(s, parallel_chains = 4, seed = 1, adapt_delta = 0.85,
                         show_exceptions = FALSE, refresh = 0,
                         iter_warmup = 500, iter_sampling = 100)
    return(r)
  }

  if (num == 6) {
    if (interactive()) {
      cat("12 indis, 1 behavior\n")
    }

    x <- generate_data(n_ids = 12, n_beh = 1,
                       beh_intercepts = c(0.5),
                       behav_types = c("count"),
                       indi_sd = 1.2, dyad_sd = 0.8)

    s <- x$standat
    r <- sociality_model(s, parallel_chains = 4, seed = 1, adapt_delta = 0.85,
                         show_exceptions = FALSE, refresh = 0,
                         iter_warmup = 500, iter_sampling = 100)
    return(r)
  }

}
