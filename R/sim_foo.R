#' run a single simulation on artificial data
#'
#' generate data, fit the model and store results to file
#'
#' @param ... arguments for \code{\link{generate_data}}
#' @param sdatalist optional list with results from  \code{\link{generate_data}}.
#'          If this is supplied, no new data will be generated in this
#'          function
#' @param indi_sd numeric, SD for gregariousness (defaults to
#'                \code{runif(1, 0.1, 3)})
#' @param dyad_sd numeric, SD for dyadic affinity (defaults to
#'                \code{runif(1, 0.1, 3)})
#' @param intercepts numeric of length two, determines the range for
#'                        intercepts. Defaults to \code{c(-3, 1)}. How many
#'                        intercepts are needed is determined inside the
#'                        function.
#' @param indi_cor,dyad_cor numeric value for correlation parameter
#'          (only works if \code{n_beh = 2}). Default is \code{NULL}, i.e.,
#'          no correlations.
#' @param n_ids see \code{\link{generate_data}}
#' @param n_beh see \code{\link{generate_data}}
#' @param behav_types see \code{\link{generate_data}}
#' @param disp_pars_gamma see \code{\link{generate_data}}
#' @param disp_pars_beta see \code{\link{generate_data}}
#' @param outpath character, output path (defaults to R's \code{\link{tempdir}})
#' @param outfilename character, output file name (defaults to a random string)
#' @param silent logical, suppress all textual output
#' @param diagnostics logical, report problems with sampling
#' @param parallel_chains number of parallel processes
#' @param chains number of total processes
#' @param refresh refresh argument
#' @param iter_warmup iter_warmup argument
#' @param iter_sampling iter_sampling argument
#' @param adapt_delta adapt_delta argument
#' @param max_treedepth max_treedepth argument
#' @param testing_mode logical, default is \code{FALSE}. If \code{TRUE}: runs
#'          the model as prior simulation, which illustrates how the function
#'          works while taking less time.
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
#' @export
#' @examples
#' \dontrun{
#' # path to results
#' # outpath <- "~/Desktop/bamosoeval" # probably doesn't exist on your computer
#' outpath <- tempdir() # default in function
#'
#' # use testing mode (faster, but only as demo)
#' for (i in 1:5) {
#'   sim_foo(silent = TRUE, outpath = outpath, testing_mode = TRUE)
#' }
#'
#' # analyse (read inidivual sim runs and collate to single data frame)
#' l <- list.files(outpath, pattern = ".rds", full.names = TRUE)
#' xdata <- sapply(l, function(x) {
#'   zz <- readRDS(x)
#'   zz$results
#' }, simplify = FALSE)
#' xdata <- do.call("rbind", xdata)
#'
#' # plot results
#' plot(xdata$indi_sd, xdata$indi_sd_est,
#'      xlab = 'true gregariousness SD', ylab = 'estimated SD')
#' abline(0, 1)
#' plot(xdata$intercept1, xdata$intercept1_est,
#'      xlab = 'true intercept (behavior 1)', ylab = 'estimated intercept (behavior 1)')
#' abline(0, 1)
#' }

# indi_sd = runif(1, 0.1, 2)
# dyad_sd = runif(1, 0.1, 2)
# n_beh = sample(3, 1)
# indi_sd = runif(2, 0.1, 2)
# dyad_sd = runif(2, 0.1, 2)
# n_beh = 2
# indi_cor <- runif(1, -0.8, 0.8)
# dyad_cor <- runif(1, -0.8, 0.8)
#
# (behav_types <- sample(c("count", "prop", "dur_gamma", "dur_beta"), n_beh, replace = TRUE))
# # (behav_types <- c("dur_gamma", "dur_beta", "dur_gamma"))
# intercepts = runif(n_beh, -2, 1)
# n_ids = sample(5:20, 1)
# disp_pars_gamma = runif(n_beh, 1, 20)
# disp_pars_beta = runif(n_beh, 5, 30)
# outpath = NULL
# outfilename = NULL
# silent = FALSE
# diagnostics = TRUE
# parallel_chains = 4; chains =4
# refresh = 200
# iter_warmup = 500; iter_sampling = 100;adapt_delta=0.9;max_treedepth=12
# testing_mode = FALSE

sim_foo <- function(...,
                    sdatalist = NULL,
                    n_ids = sample(5:60, 1),
                    indi_sd = runif(1, 0.1, 3),
                    dyad_sd = runif(1, 0.1, 3),
                    n_beh = sample(3, 1),
                    indi_cor = NULL,
                    dyad_cor = NULL,
                    behav_types = NULL,
                    intercepts = runif(n_beh, -3, 1),
                    disp_pars_gamma = runif(n_beh, 1, 20),
                    disp_pars_beta = runif(n_beh, 5, 30),
                    outpath = NULL,
                    outfilename = NULL,
                    silent = TRUE,
                    diagnostics = TRUE,
                    chains = 4,
                    parallel_chains = 1,
                    refresh = 0,
                    iter_warmup = 2000,
                    iter_sampling = 1000,
                    adapt_delta = 0.98,
                    max_treedepth = 14,
                    testing_mode = FALSE
                    ) {

  if (testing_mode) {
    if (interactive()) {
      message("you are in testing mode: the posteriors are not conditioned on the data")
    }
  }

  if (is.null(sdatalist)) {
    if (is.null(behav_types)) {
      behav_types <- sample(c("count", "prop", "dur_gamma", "dur_beta"), n_beh, replace = TRUE)
    }

    if (!is.null(indi_cor)) {
      if (n_beh != 2) stop("correlations for gregariousness currently requires exactly n_beh=2")
      if (length(indi_sd) != 2) stop("when using correlations for gregariousness currently requires vector of length 2 for indi_sd")
      indi_sd <- matrix(c(indi_sd[1], indi_cor, indi_cor, indi_sd[2]), ncol = 2)
    }
    if (!is.null(dyad_cor)) {
      if (n_beh != 2) stop("correlations for affinity currently requires exactly n_beh=2")
      if (length(dyad_sd) != 2) stop("when using correlations for affinity currently requires vector of length 2 for dyad_sd")
      dyad_sd <- matrix(c(dyad_sd[1], dyad_cor, dyad_cor, dyad_sd[2]), ncol = 2)
    }

    # occasionally, an empty matrix is produced (mostly in the gamma setting)
    good_to_go <- FALSE
    while (!good_to_go) {
      sdata <- generate_data(n_beh = n_beh,
                             n_ids = n_ids,
                             behav_types = behav_types,
                             beh_intercepts = intercepts,
                             indi_sd = indi_sd,
                             dyad_sd = dyad_sd,
                             disp_pars_gamma = disp_pars_gamma,
                             disp_pars_beta = disp_pars_beta,
                             # ...
      )
      if (!sdata$input_data$empty) good_to_go <- TRUE
    }

    # ignore data when matrices with only 0s were generated
    # this check should be redundant, but better safe than sorry
    if (sdata$input_data$empty) {
      cat("empty matrix generated: skipped\n")
      return(NULL)
    }
  } else {
    sdata <- sdatalist$sdata
    behav_types <- names(sdata$standat$behav_types)[1:sdata$standat$n_beh]
    indi_cor <- NULL
    if (isTRUE(ncol(sdata$input_data$indi_sd) == 2)) {
      indi_cor <- sdata$input_data$indi_sd[1, 2]
    }
    dyad_cor <- NULL
    if (isTRUE(ncol(sdata$input_data$dyad_sd) == 2)) {
      dyad_cor <- sdata$input_data$dyad_sd[1, 2]
    }
    intercepts <- sdata$input_data$beh_intercepts
  }


  # fit the model
  seed <- sample(1e+7, 1)

  if (silent) {
    r <- sociality_model(standat = sdata$standat,
                         silent = TRUE,
                         show_messages = FALSE,
                         diagnostics = NULL,
                         parallel_chains = parallel_chains,
                         refresh = refresh,
                         iter_warmup = iter_warmup,
                         iter_sampling = iter_sampling,
                         adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth,
                         seed = seed,
                         prior_sim = testing_mode,
                         chains = chains
                         )
  } else {
    r <- sociality_model(standat = sdata$standat,
                         silent = FALSE,
                         show_messages = FALSE,
                         diagnostics = NULL,
                         parallel_chains = parallel_chains,
                         refresh = refresh,
                         iter_warmup = iter_warmup,
                         iter_sampling = iter_sampling,
                         adapt_delta = adapt_delta,
                         max_treedepth = max_treedepth,
                         seed = seed,
                         prior_sim = testing_mode,
                         chains = chains
                         )
  }


  # extract results
  xres_draws <- r$mod_res$draws(format = "draws_matrix")
  # clean up draws (remove irrelevant pars) (for simulations anyway)
  xres_draws <- xres_draws[, !grepl("log_lik", colnames(xres_draws))]
  xres_draws <- xres_draws[, !grepl("interactions_pred_cont", colnames(xres_draws))]
  xres_draws <- xres_draws[, !grepl("interactions_pred", colnames(xres_draws))]
  xres_draws <- xres_draws[, !grepl("scaled_indi_sums", colnames(xres_draws))]
  xres_draws <- xres_draws[, !grepl("dyad_soc_vals_z", colnames(xres_draws))]
  xres_draws <- xres_draws[, !grepl("indi_soc_vals_z", colnames(xres_draws))]
  colnames(xres_draws)



  if (!"beh_intercepts[1]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "beh_intercepts[1]" = NA)
  }
  if (!"beh_intercepts[2]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "beh_intercepts[2]" = NA)
  }
  if (!"beh_intercepts[3]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "beh_intercepts[3]" = NA)
  }

  if (!"shapes_gamma[1]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "shapes_gamma[1]" = NA)
  }
  if (!"shapes_gamma[2]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "shapes_gamma[2]" = NA)
  }
  if (!"shapes_gamma[3]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "shapes_gamma[3]" = NA)
  }

  if (!"shapes_beta[1]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "shapes_beta[1]" = NA)
  }
  if (!"shapes_beta[2]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "shapes_beta[2]" = NA)
  }
  if (!"shapes_beta[3]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "shapes_beta[3]" = NA)
  }

  # move draws to correct position for shape parameters
  # sdata$standat$beta_shape_n
  # sdata$standat$beta_shape_pos
  #
  # sdata$standat$gamma_shape_n
  # sdata$standat$gamma_shape_pos

  checksum1 <- sum(round(xres_draws[, grepl("shapes_gamma", colnames(xres_draws))]), na.rm = TRUE)
  if (isTRUE(sdata$standat$gamma_shape_pos[3] == 1)) {
    xres_draws[, "shapes_gamma[3]"] <- xres_draws[, "shapes_gamma[1]"]
    xres_draws[, "shapes_gamma[1]"] <- NA
  }
  if (isTRUE(sdata$standat$gamma_shape_pos[3] == 2)) {
    xres_draws[, "shapes_gamma[3]"] <- xres_draws[, "shapes_gamma[2]"]
    xres_draws[, "shapes_gamma[2]"] <- NA
  }
  if (isTRUE(sdata$standat$gamma_shape_pos[2] == 1)) {
    xres_draws[, "shapes_gamma[2]"] <- xres_draws[, "shapes_gamma[1]"]
    xres_draws[, "shapes_gamma[1]"] <- NA
  }
  checksum2 <- sum(round(xres_draws[, grepl("shapes_gamma", colnames(xres_draws))]), na.rm = TRUE)

  checksum3 <- sum(round(xres_draws[, grepl("shapes_beta", colnames(xres_draws))]), na.rm = TRUE)
  if (isTRUE(sdata$standat$beta_shape_pos[3] == 1)) {
    xres_draws[, "shapes_beta[3]"] <- xres_draws[, "shapes_beta[1]"]
    xres_draws[, "shapes_beta[1]"] <- NA
  }
  if (isTRUE(sdata$standat$beta_shape_pos[3] == 2)) {
    xres_draws[, "shapes_beta[3]"] <- xres_draws[, "shapes_beta[2]"]
    xres_draws[, "shapes_beta[2]"] <- NA
  }
  if (isTRUE(sdata$standat$beta_shape_pos[2] == 1)) {
    xres_draws[, "shapes_beta[2]"] <- xres_draws[, "shapes_beta[1]"]
    xres_draws[, "shapes_beta[1]"] <- NA
  }
  checksum4 <- sum(round(xres_draws[, grepl("shapes_beta", colnames(xres_draws))]), na.rm = TRUE)

  msg <- "rearranging draws error"
  if (isTRUE(behav_types[1] == "dur_gamma")) if (any(is.na(xres_draws[, "shapes_gamma[1]"]))) stop(msg)
  if (isTRUE(behav_types[2] == "dur_gamma")) if (any(is.na(xres_draws[, "shapes_gamma[2]"]))) stop(msg)
  if (isTRUE(behav_types[3] == "dur_gamma")) if (any(is.na(xres_draws[, "shapes_gamma[3]"]))) stop(msg)
  if (isTRUE(behav_types[1] == "dur_beta")) if (any(is.na(xres_draws[, "shapes_beta[1]"]))) stop(msg)
  if (isTRUE(behav_types[2] == "dur_beta")) if (any(is.na(xres_draws[, "shapes_beta[2]"]))) stop(msg)
  if (isTRUE(behav_types[3] == "dur_beta")) if (any(is.na(xres_draws[, "shapes_beta[3]"]))) stop(msg)

  # unify non-cor with cor models
  if ("indi_soc_sd" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "indi_soc_sd[1]" = xres_draws[, "indi_soc_sd"], "indi_soc_sd[2]" = NA)
    xres_draws <- xres_draws[, -c(which(colnames(xres_draws) == "indi_soc_sd"))]
  }
  if ("dyad_soc_sd" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "dyad_soc_sd[1]" = xres_draws[, "dyad_soc_sd"], "dyad_soc_sd[2]" = NA)
    xres_draws <- xres_draws[, -c(which(colnames(xres_draws) == "dyad_soc_sd"))]
  }
  if (!"cors_indi[1]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "cors_indi[1]" = NA)
  }
  if (!"cors_dyad[1]" %in% colnames(xres_draws)) {
    xres_draws <- cbind(xres_draws, "cors_dyad[1]" = NA)
  }



  # quantiles
  probs <- c(0.1, 0.9)
  qres80 <- data.frame(
    bi1l = quantile(c(xres_draws[, "beh_intercepts[1]"]), probs = probs[1], na.rm = TRUE),
    bi1u = quantile(c(xres_draws[, "beh_intercepts[1]"]), probs = probs[2], na.rm = TRUE),
    bi2l = quantile(c(xres_draws[, "beh_intercepts[2]"]), probs = probs[1], na.rm = TRUE),
    bi2u = quantile(c(xres_draws[, "beh_intercepts[2]"]), probs = probs[2], na.rm = TRUE),
    bi3l = quantile(c(xres_draws[, "beh_intercepts[3]"]), probs = probs[1], na.rm = TRUE),
    bi3u = quantile(c(xres_draws[, "beh_intercepts[3]"]), probs = probs[2], na.rm = TRUE),

    sg1l = quantile(c(xres_draws[, "shapes_gamma[1]"]), probs = probs[1], na.rm = TRUE),
    sg1u = quantile(c(xres_draws[, "shapes_gamma[1]"]), probs = probs[2], na.rm = TRUE),
    sg2l = quantile(c(xres_draws[, "shapes_gamma[2]"]), probs = probs[1], na.rm = TRUE),
    sg2u = quantile(c(xres_draws[, "shapes_gamma[2]"]), probs = probs[2], na.rm = TRUE),
    sg3l = quantile(c(xres_draws[, "shapes_gamma[3]"]), probs = probs[1], na.rm = TRUE),
    sg3u = quantile(c(xres_draws[, "shapes_gamma[3]"]), probs = probs[2], na.rm = TRUE),

    sb1l = quantile(c(xres_draws[, "shapes_beta[1]"]), probs = probs[1], na.rm = TRUE),
    sb1u = quantile(c(xres_draws[, "shapes_beta[1]"]), probs = probs[2], na.rm = TRUE),
    sb2l = quantile(c(xres_draws[, "shapes_beta[2]"]), probs = probs[1], na.rm = TRUE),
    sb2u = quantile(c(xres_draws[, "shapes_beta[2]"]), probs = probs[2], na.rm = TRUE),
    sb3l = quantile(c(xres_draws[, "shapes_beta[3]"]), probs = probs[1], na.rm = TRUE),
    sb3u = quantile(c(xres_draws[, "shapes_beta[3]"]), probs = probs[2], na.rm = TRUE),

    indil = quantile(c(xres_draws[, "indi_soc_sd[1]"]), probs = probs[1]),
    indiu = quantile(c(xres_draws[, "indi_soc_sd[1]"]), probs = probs[2]),
    dyadl = quantile(c(xres_draws[, "dyad_soc_sd[1]"]), probs = probs[1]),
    dyadu = quantile(c(xres_draws[, "dyad_soc_sd[1]"]), probs = probs[2]),

    indil2 = quantile(c(xres_draws[, "indi_soc_sd[2]"]), probs = probs[1], na.rm = TRUE),
    indiu2 = quantile(c(xres_draws[, "indi_soc_sd[2]"]), probs = probs[2], na.rm = TRUE),
    dyadl2 = quantile(c(xres_draws[, "dyad_soc_sd[2]"]), probs = probs[1], na.rm = TRUE),
    dyadu2 = quantile(c(xres_draws[, "dyad_soc_sd[2]"]), probs = probs[2], na.rm = TRUE),

    indicorl = quantile(c(xres_draws[, "cors_indi[1]"]), probs = probs[1], na.rm = TRUE),
    indicoru = quantile(c(xres_draws[, "cors_indi[1]"]), probs = probs[2], na.rm = TRUE),
    dyadcorl = quantile(c(xres_draws[, "cors_dyad[1]"]), probs = probs[1], na.rm = TRUE),
    dyadcoru = quantile(c(xres_draws[, "cors_dyad[1]"]), probs = probs[2], na.rm = TRUE)
  )

  probs <- c(0.055, 0.945)
  qres89 <- data.frame(
    bi1l = quantile(c(xres_draws[, "beh_intercepts[1]"]), probs = probs[1], na.rm = TRUE),
    bi1u = quantile(c(xres_draws[, "beh_intercepts[1]"]), probs = probs[2], na.rm = TRUE),
    bi2l = quantile(c(xres_draws[, "beh_intercepts[2]"]), probs = probs[1], na.rm = TRUE),
    bi2u = quantile(c(xres_draws[, "beh_intercepts[2]"]), probs = probs[2], na.rm = TRUE),
    bi3l = quantile(c(xres_draws[, "beh_intercepts[3]"]), probs = probs[1], na.rm = TRUE),
    bi3u = quantile(c(xres_draws[, "beh_intercepts[3]"]), probs = probs[2], na.rm = TRUE),

    sg1l = quantile(c(xres_draws[, "shapes_gamma[1]"]), probs = probs[1], na.rm = TRUE),
    sg1u = quantile(c(xres_draws[, "shapes_gamma[1]"]), probs = probs[2], na.rm = TRUE),
    sg2l = quantile(c(xres_draws[, "shapes_gamma[2]"]), probs = probs[1], na.rm = TRUE),
    sg2u = quantile(c(xres_draws[, "shapes_gamma[2]"]), probs = probs[2], na.rm = TRUE),
    sg3l = quantile(c(xres_draws[, "shapes_gamma[3]"]), probs = probs[1], na.rm = TRUE),
    sg3u = quantile(c(xres_draws[, "shapes_gamma[3]"]), probs = probs[2], na.rm = TRUE),

    sb1l = quantile(c(xres_draws[, "shapes_beta[1]"]), probs = probs[1], na.rm = TRUE),
    sb1u = quantile(c(xres_draws[, "shapes_beta[1]"]), probs = probs[2], na.rm = TRUE),
    sb2l = quantile(c(xres_draws[, "shapes_beta[2]"]), probs = probs[1], na.rm = TRUE),
    sb2u = quantile(c(xres_draws[, "shapes_beta[2]"]), probs = probs[2], na.rm = TRUE),
    sb3l = quantile(c(xres_draws[, "shapes_beta[3]"]), probs = probs[1], na.rm = TRUE),
    sb3u = quantile(c(xres_draws[, "shapes_beta[3]"]), probs = probs[2], na.rm = TRUE),

    indil = quantile(c(xres_draws[, "indi_soc_sd[1]"]), probs = probs[1]),
    indiu = quantile(c(xres_draws[, "indi_soc_sd[1]"]), probs = probs[2]),
    dyadl = quantile(c(xres_draws[, "dyad_soc_sd[1]"]), probs = probs[1]),
    dyadu = quantile(c(xres_draws[, "dyad_soc_sd[1]"]), probs = probs[2]),

    indil2 = quantile(c(xres_draws[, "indi_soc_sd[2]"]), probs = probs[1], na.rm = TRUE),
    indiu2 = quantile(c(xres_draws[, "indi_soc_sd[2]"]), probs = probs[2], na.rm = TRUE),
    dyadl2 = quantile(c(xres_draws[, "dyad_soc_sd[2]"]), probs = probs[1], na.rm = TRUE),
    dyadu2 = quantile(c(xres_draws[, "dyad_soc_sd[2]"]), probs = probs[2], na.rm = TRUE),

    indicorl = quantile(c(xres_draws[, "cors_indi[1]"]), probs = probs[1], na.rm = TRUE),
    indicoru = quantile(c(xres_draws[, "cors_indi[1]"]), probs = probs[2], na.rm = TRUE),
    dyadcorl = quantile(c(xres_draws[, "cors_dyad[1]"]), probs = probs[1], na.rm = TRUE),
    dyadcoru = quantile(c(xres_draws[, "cors_dyad[1]"]), probs = probs[2], na.rm = TRUE)

  )

  cor_indi <- NA
  cor_dyad <- NA
  if (sdata$standat$n_cors == 0) {
    indi_columns <- grepl("indi_soc_vals[", colnames(xres_draws), fixed = TRUE)
    cor_indi <- cor(sdata$input_data$indi_data$indi_soc_vals,
                    apply(xres_draws[, indi_columns], 2, median))
    dyad_columns <- grepl("dyad_soc_vals[", colnames(xres_draws), fixed = TRUE)
    cor_dyad <- cor(sdata$input_data$dyad_data$dyad_soc_vals,
                    apply(xres_draws[, dyad_columns], 2, median))
  }

  d <- r$mod_res$diagnostic_summary(quiet = TRUE)



  if (is.null(outfilename)) {
    fn <- paste(sample(c(letters, 0:9), 10, TRUE), collapse = "")
    fn <- paste0(fn, ".rds")
    if (!is.null(sdata$input_data$filename)) {
      fn <- sdata$input_data$filename
    }
  }

  # collect relevant info
  ## width of posteriors for soc components
  # w1 <- quantile(c(xres_draws[, "indi_soc_sd"]), 0.0945) - quantile(c(xres_draws[, "indi_soc_sd"]), 0.055)
  # w2 <- quantile(c(xres_draws[, "dyad_soc_sd"]), 0.0945) - quantile(c(xres_draws[, "dyad_soc_sd"]), 0.055)

  n_beh <- sdata$standat$n_beh
  xres_results <- data.frame(label = fn,
                             n_ids = sdata$standat$n_ids,
                             n_beh = sdata$standat$n_beh,
                             indi_sd1 = sdata$input_data$indi_sd[1],
                             indi_sd1_est = median(xres_draws[, "indi_soc_sd[1]"]),
                             indi_sd2 = sdata$input_data$indi_sd[4], # 4 in a 2x2 matrix, uhh...
                             indi_sd2_est = median(xres_draws[, "indi_soc_sd[2]"]),
                             # width_indisd = w1,
                             dyad_sd1 = sdata$input_data$dyad_sd[1],
                             dyad_sd1_est = median(xres_draws[, "dyad_soc_sd[1]"]),
                             dyad_sd2 = sdata$input_data$dyad_sd[4],
                             dyad_sd2_est = median(xres_draws[, "dyad_soc_sd[2]"]),
                             # width_dyadsd = w2,
                             indi_cor = ifelse(is.null(indi_cor), NA, indi_cor),
                             indi_cor_est = median(xres_draws[, "cors_indi[1]"]),
                             dyad_cor = ifelse(is.null(dyad_cor), NA, dyad_cor),
                             dyad_cor_est = median(xres_draws[, "cors_dyad[1]"]),

                             intercept1 = intercepts[1],
                             intercept1_est = median(xres_draws[, "beh_intercepts[1]"]),
                             intercept2 = intercepts[2],
                             intercept2_est = ifelse(n_beh >= 2, median(xres_draws[, "beh_intercepts[2]"]), NA),
                             intercept3 = intercepts[3],
                             intercept3_est = ifelse(n_beh >= 3, median(xres_draws[, "beh_intercepts[3]"]), NA),
                             shapegamma1 = NA,
                             shapegamma1_est = NA,
                             shapegamma2 = NA,
                             shapegamma2_est = NA,
                             shapegamma3 = NA,
                             shapegamma3_est = NA,
                             shapebeta1 = NA,
                             shapebeta1_est = NA,
                             shapebeta2 = NA,
                             shapebeta2_est = NA,
                             shapebeta3 = NA,
                             shapebeta3_est = NA,
                             blups_cor_indi = cor_indi,
                             blups_cor_dyad = cor_dyad,
                             tot_dur = round(r$mod_res$time()$total, 2),
                             n_samples = nrow(xres_draws),
                             divergent = sum(d$num_divergent),
                             treedepth = sum(d$num_max_treedepth),
                             ebfmi = sum(d$ebfmi < 0.2),
                             seed = seed,
                             behav_types = paste(sort(behav_types), collapse = ";"),
                             b1 = sdata$input_data$behav_type[1],
                             b2 = sdata$input_data$behav_type[2],
                             b3 = sdata$input_data$behav_type[3]

  )

  # shape parameters
  if (sdata$standat$gamma_shape_n > 0) {
    cnt <- 1
    for (i in seq_len(sdata$standat$gamma_shape_n)) {
      pos <- which(sdata$standat$gamma_shape_pos == cnt)
      aux <- c(xres_draws[, paste0("shapes_gamma[", pos, "]")])
      xres_results[1, paste0("shapegamma", pos, "_est")] <- median(aux)
      xres_results[1, paste0("shapegamma", pos)] <- sdata$input_data$disp_pars_gamma[pos]
      cnt <- cnt + 1
    }
  }
  if (sdata$standat$beta_shape_n > 0) {
    cnt <- 1
    for (i in seq_len(sdata$standat$beta_shape_n)) {
      pos <- which(sdata$standat$beta_shape_pos == cnt)
      aux <- c(xres_draws[, paste0("shapes_beta[", pos, "]")])
      xres_results[1, paste0("shapebeta", pos, "_est")] <- median(aux)
      xres_results[1, paste0("shapebeta", pos)] <- sdata$input_data$disp_pars_beta[pos]
      cnt <- cnt + 1
    }
  }

  if (!is.null(sdata$input_data$filename)) {
    fn <- sdata$input_data$filename
  } else {
    if (is.null(outfilename)) {
      fn <- paste(sample(c(letters, 0:9), 10, TRUE), collapse = "")
      fn <- paste0(fn, ".rds")
    }
  }


  if (is.null(outpath)) {
    op <- normalizePath(tempdir())
  } else {
    op <- normalizePath(outpath)
  }

  outloc <- file.path(op, fn)
  res <- list(sdata = sdata, # xres = xres_summary,
              results = xres_results,
              qres80 = qres80, # quantiles for key pars
              qres89 = qres89, # quantiles for key pars
              diagnostics = d,
              draws = xres_draws,
              seed = seed)
  suppressWarnings(saveRDS(res, file = outloc))

  dres <- "all good"
  if (any(d$num_divergent > 0)) dres <- "sampling problems (check '$diagnostics')"
  if (any(d$num_max_treedepth > 0)) dres <- "sampling problems (check '$diagnostics')"
  if (any(d$ebfmi < 0.2)) dres <- "sampling problems (check '$diagnostics')"
  if (!silent) {
    if (diagnostics) cat(dres, "\n")
  }

  # return NULL or file path?
  if (silent) {
    return(invisible(outloc))
  }
  outloc
}


#' create simulated data sets
#'
#' @inheritParams sim_foo
#'
#' @returns a list
#' @export
#'

sim_datasets <- function(n_ids = sample(5:60, 1),
                         indi_sd = runif(1, 0.1, 3),
                         dyad_sd = runif(1, 0.1, 3),
                         n_beh = sample(3, 1),
                         indi_cor = NULL,
                         dyad_cor = NULL,
                         behav_types = NULL,
                         intercepts = runif(n_beh, -3, 1),
                         disp_pars_gamma = runif(n_beh, 1, 20),
                         disp_pars_beta = runif(n_beh, 5, 30),
                         outpath = NULL,
                         outfilename = NULL,
                         silent = TRUE
) {

  if (is.null(behav_types)) {
    behav_types <- sample(c("count", "prop", "dur_gamma", "dur_beta"), n_beh, replace = TRUE)
  }

  if (!is.null(indi_cor)) {
    if (n_beh != 2) stop("correlations for gregariousness currently requires exactly n_beh=2")
    if (length(indi_sd) != 2) stop("when using correlations for gregariousness currently requires vector of length 2 for indi_sd")
    indi_sd <- matrix(c(indi_sd[1], indi_cor, indi_cor, indi_sd[2]), ncol = 2)
  }
  if (!is.null(dyad_cor)) {
    if (n_beh != 2) stop("correlations for affinity currently requires exactly n_beh=2")
    if (length(dyad_sd) != 2) stop("when using correlations for affinity currently requires vector of length 2 for dyad_sd")
    dyad_sd <- matrix(c(dyad_sd[1], dyad_cor, dyad_cor, dyad_sd[2]), ncol = 2)
  }

  # occasionally, an empty matrix is produced (mostly in the gamma setting)
  good_to_go <- FALSE
  while (!good_to_go) {
    sdata <- generate_data(n_beh = n_beh,
                           n_ids = n_ids,
                           behav_types = behav_types,
                           beh_intercepts = intercepts,
                           indi_sd = indi_sd,
                           dyad_sd = dyad_sd,
                           disp_pars_gamma = disp_pars_gamma,
                           disp_pars_beta = disp_pars_beta,
                           # ...
    )
    if (!sdata$input_data$empty) good_to_go <- TRUE
  }

  # ignore data when matrices with only 0s were generated
  # this check should be redundant, but better safe than sorry
  if (sdata$input_data$empty) {
    cat("empty matrix generated: skipped\n")
    return(NULL)
  }

  if (is.null(outfilename)) {
    fn <- paste(sample(c(letters, 0:9), 10, TRUE), collapse = "")
    fn <- paste0(fn, ".rds")
  } else {
    fn <- outfilename
  }

  if (is.null(outpath)) {
    op <- normalizePath(tempdir())
  } else {
    op <- normalizePath(outpath)
  }

  sdata$input_data$filename <- fn

  outloc <- file.path(op, fn)
  res <- list(sdata = sdata)
  if (file.exists(outloc)) {
    message("skipped ", shQuote(fn), " because it already exists")
  } else {
    suppressWarnings(saveRDS(res, file = outloc))
  }


  if (silent) {
    return(invisible(outloc))
  }
  outloc
}
