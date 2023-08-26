#' run a single simulation on artificial data
#'
#' fit the model with simulated data and store results to file
#'
#' @param ... arguments for \code{\link{generate_data}}
#' @param indi_sd numeric, SD for gregariousness (defaults to
#'                \code{runif(1, 0.1, 3)})
#' @param dyad_sd numeric, SD for dyadic affinity (defaults to
#'                \code{runif(1, 0.1, 3)})
#' @param intercept_range numeric of length two, determines the range for
#'                        intercepts. Defaults to \code{c(-3, 0)}. How many
#'                        intercepts are needed is determined inside the
#'                        function.
#' @param outpath character, output path (defaults to R's \code{\link{tempdir}})
#' @param outfilename character, output file name (defaults to a random string)
#' @param silent logical, suppress all textual output
#' @param diagnostics logical, report problems with sampling
#' @param parallel_chains number of parallel processes
#' @param refresh refresh argument
#' @param iter_warmup iter_warmup argument
#' @param iter_sampling iter_sampling argument
#'
#' @importFrom utils capture.output
#'
#' @return writes an rds file to the specified location, which contains a list
#'         with the data generation settings, generated data,
#'         fitted model's draws, and seed used during the fitting.
#'         Most importantly, the list contains a one-line data frame
#'         (\code{$results}) with the key results used for downstream
#'         analyses (see Details section).
#'
#' @details
#' The function starts by generating interaction data across variable settings
#' (see \code{\link{generate_data}}). Then it fits the model to these data.
#' From the model results it extracts the following quantities, which are in
#' the \code{$results} list item.
#'
#' * number of individuals
#' * number of behaviors
#' * true and estimated (posterior median) gregariousness SD
#' * true and estimated (posterior median) affinity SD
#' * correlation between true and estimated (posterior medians) of individual
#' gregariousness values
#' * correlation between true and estimated (posterior medians) of dyadic
#' affinity values
#' * number of posterior samples
#' * sampler diagnostics (divergent transitions etc)
#'
#'
#'
#' @export
#'


sim_foo <- function(...,
                    indi_sd = runif(1, 0.1, 3),
                    dyad_sd = runif(1, 0.1, 3),
                    intercept_range = c(-3, 0),
                    outpath = NULL,
                    outfilename = NULL,
                    silent = TRUE,
                    diagnostics = TRUE,
                    parallel_chains = 1,
                    refresh = 0,
                    iter_warmup = 2000,
                    iter_sampling = 1000) {

  xx <- capture.output(mod <- suppressMessages(get_model()))
  rm(xx)

  n_beh <- sample(1:3, 1)
  behav_types <- sample(c("count", "prop", "dur_gamma", "dur_beta"), n_beh, replace = TRUE)
  intercepts <- runif(n_beh, intercept_range[1], intercept_range[2])

  # occasionally, an empty matrix is produced (mostly in the gamma setting)
  good_to_go <- FALSE
  while (!good_to_go) {
    sdata <- generate_data(n_beh = n_beh,
                           behav_types = behav_types,
                           beh_intercepts = intercepts,
                           indi_sd = indi_sd,
                           dyad_sd = dyad_sd,
                           ...)
    if (!sdata$input_data$empty) good_to_go <- TRUE
  }


  # ignore data when matrices with only 0s were generated
  # this check should be redundant, but better safe than sorry
  if (sdata$input_data$empty) {
    cat("empty matrix generated: skipped\n")
    return(NULL)
  }

  # fit the model
  seed <- sample(1e+7, 1)
  xx <- capture.output(
    xres <- suppressMessages(mod$sample(data = sdata$standat,
                                        show_messages = FALSE,
                                        diagnostics = NULL,
                                        parallel_chains = parallel_chains,
                                        refresh = refresh,
                                        iter_warmup = iter_warmup,
                                        iter_sampling = iter_sampling,
                                        adapt_delta = 0.98,
                                        max_treedepth = 14,
                                        seed = seed))
  )

  # extract results
  # not really needed: xres_summary <- xres$summary()
  xres_draws <- xres$draws(format = "draws_matrix")
  indi_columns <- grepl("indi_soc_vals[", colnames(xres_draws), fixed = TRUE)
  cor_indi <- cor(sdata$input_data$indi_soc_vals,
                  apply(xres_draws[, indi_columns], 2, median))
  dyad_columns <- grepl("dyad_soc_vals[", colnames(xres_draws), fixed = TRUE)
  cor_dyad <- cor(sdata$input_data$dyad_soc_vals,
                  apply(xres_draws[, dyad_columns], 2, median))

  d <- xres$diagnostic_summary(quiet = TRUE)

  # collect relevant info
  xres_results <- data.frame(n_ids = sdata$standat$n_ids,
                             n_beh = n_beh,
                             indi_sd = sdata$input_data$indi_sd,
                             indi_sd_est = median(xres_draws[, "indi_soc_sd"]),
                             dyad_sd = sdata$input_data$dyad_sd,
                             dyad_sd_est = median(xres_draws[, "dyad_soc_sd"]),
                             cor_indi = cor_indi,
                             cor_dyad = cor_dyad,
                             tot_dur = round(xres$time()$total, 2),
                             n_samples = nrow(xres_draws),
                             divergent = sum(d$num_divergent),
                             treedepth = sum(d$num_max_treedepth),
                             ebfmi = sum(d$ebfmi < 0.2)
  )

  if (is.null(outfilename)) {
    fn <- paste(sample(c(letters, 0:9), 10, TRUE), collapse = "")
    fn <- paste0(fn, ".rds")
  }
  if (is.null(outpath)) {
    op <- normalizePath(tempdir())
  } else {
    op <- normalizePath(outpath)
  }
  outloc <- file.path(op, fn)
  res <- list(sdata = sdata, # xres = xres_summary,
              results = xres_results,
              diagnostics = d,
              draws = xres_draws,
              seed = seed)
  suppressWarnings(saveRDS(res, file = outloc))

  dres <- "all good"
  if (any(d$num_divergent > 0)) dres <- "sampling problems"
  if (any(d$num_max_treedepth > 0)) dres <- "sampling problems"
  if (any(d$ebfmi < 0.2)) dres <- "sampling problems"
  if (!silent) {
    if (diagnostics) cat(dres, "\n")
  }

  # return NULL or file path?
  if (silent) {
    return(invisible(outloc))
  }
  outloc
}
