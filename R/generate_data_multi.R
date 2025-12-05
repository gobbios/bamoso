#' generate multiple groups with features from a underlying mechanisms
#'
#' please consider this still EXPERIMENTAL!
#'
#' @param groups list with group compositions (see \code{\link{generate_groups}})
#' @param indisds,dyadsds vector with SD
#' @param bintercepts vector with intercepts
#'
#' @returns a list with matrices that can be passed to
#' \code{\link{make_stan_data_from_matrices_multi}}
#'
#' @details
#' The function is a wrapper for \code{\link{generat_data}}.
#'
#'
#' @examples
#' n <- 3 # number of groups to generate
#' # just assign group size and id codes
#' groups <- generate_groups(n_groups = n, min_nids = 4, max_nids = 5)
#' # set true parameter value
#' (indisds <- round(exp(rnorm(n, -1, 0.8)), 1))
#' (dyadsds <- round(exp(rnorm(n, -0.5, 0.6)), 1))
#' (bintercepts <- round(rnorm(n, -2, 1), 1))
#' generate_data_multi(groups, indisds, dyadsds, bintercepts)

# n <- 5
# groups <- generate_groups(n_groups = n, min_nids = 5, max_nids = 12)
# (indisds <- round(exp(rnorm(n, -1, 0.8)), 1))
# (dyadsds <- round(exp(rnorm(n, -0.5, 0.6)), 1))
# (bintercepts <- round(rnorm(n, -2, 1), 1))

generate_data_multi <- function(groups, indisds, dyadsds, bintercepts) {
  blist <- vector("list", length(groups))
  olist <- vector("list", length(groups))
  names(blist) <- names(groups)
  names(olist) <- names(groups)

  i=1
  for (i in seq_along(groups)) {
    xx <- generate_data(n_ids = length(groups[[i]]), n_beh = 1, behav_types = "count",
                        indi_sd = indisds[i], dyad_sd = dyadsds[i],
                        beh_intercepts = bintercepts[i], count_obseff = c(0.5, 5))
    x <- xx$processed$interaction_matrices[[1]]
    colnames(x) <- rownames(x) <- groups[[i]]
    blist[[i]] <- x
    x <- xx$processed$obseff_matrices[[1]]
    colnames(x) <- rownames(x) <- groups[[i]]
    olist[[i]] <- x
  }


  list(bmats = blist, omats = olist)
}
