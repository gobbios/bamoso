#' summarize dyadicmodel object
#'
#' @param object an object of class \code{"dyadicmodel"},
#'               usually the result of a call to \code{\link{sociality_model}}
#' @param ... further arguments passed to or from other methods (ignored)
#' @importFrom stats median quantile
#' @author Christof Neumann
#'
#' @export

summary.dyadicmodel <- function(object, ...) {
  d <- object$mod_res$draws(c("indi_soc_sd", "dyad_soc_sd"),
                            format = "draws_matrix")
  qs_i <- sprintf("%.2f", quantile(d[, "indi_soc_sd"], c(0.055, 0.5, 0.945)))
  qs_d <- sprintf("%.2f", quantile(d[, "dyad_soc_sd"], c(0.055, 0.5, 0.945)))

  s <- paste0("Model of dyadic interactions from ", object$standat$n_ids,
              " individuals (", object$standat$n_dyads, " dyads)")
  cat(s, "\n")
  if ("removed_dyads" %in% names(object$standat)) {
    cat("(removed", nrow(object$standat$removed_dyads), "dyads with NA values)\n")
  }

  nchains <- length(object$mod_res$metadata()$id)
  sperchain <- object$mod_res$metadata()$iter_sampling
  cat("number of post-warmup samples:", nchains * sperchain, "\n")

  cat("--------------------\n")
  cat("sociality estimates:\n")
  cat("  individual gregariousness SD (median [89% CI]):",
      qs_i[2], paste0("[", qs_i[1], " - ", qs_i[3], "]\n"))
  cat("  pairwise affinity SD (median [89% CI]):        ",
      qs_d[2], paste0("[", qs_d[1], " - ", qs_d[3], "]\n"))
  cat("--------------------\n")

  cat(object$standat$n_beh, "behavior(s) included:\n")
  for (i in seq_len(object$standat$n_beh)) {
    s <- paste0("(", i, "): ",
                shQuote(names(object$standat$beh_names)[i]), " (type: ",
                shQuote(names(object$standat$behav_types)[i]), ")")
    cat(s, "\n")
  }

  cat("--------------------\n")
  diagnostics <- object$mod_res$diagnostic_summary()
  if (all(diagnostics$num_divergent == 0) &&
      all(diagnostics$num_max_treedepth == 0) &&
      all(diagnostics$ebfmi > 0.2)) {
    diags <- TRUE
    msg <- "no obvious sampling issues detected"
  } else {
    diags <- FALSE
    msg <- "sampling issues detected (see below)"
  }
    if (diags) cat(msg, "\n")
  if (!diags) {
    cat(msg, "\n")
    if (interactive()) {
      print(diagnostics)
    }
  }
  cat("\n")

}
