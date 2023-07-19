#' plot posteriors of individual and dyadic sociality SDs
#'
#' @param mod_res model result from \code{\link{sociality_model}}
#' @param xlim limits for horizontal axis
#' @param labs named list with labels for the individual- and
#'             dyad-specific component
#' @param do_legend logical, should the legend be plotted
#'                  (default is \code{TRUE})
#' @param do_dyadic logical, should the posterior for the dyad-level
#'                  estimate be plotted (default is \code{TRUE})
#' @param do_indi logical, should the posterior for the individual-level
#'                estimate be plotted (default is \code{TRUE})
#'
#' @return a plot
#' @export
#'
#' @importFrom brms rstudent_t
#' @importFrom grDevices hcl.colors
#' @importFrom graphics legend polygon

sociality_plot <- function(mod_res,
                           xlim = c(0, 4),
                           labs = list(indi = "gregariousness",
                                       dyad = "affinity"),
                           do_legend = TRUE,
                           do_dyadic = TRUE,
                           do_indi = TRUE) {

  mod_res <- mod_res$mod_res

  xcols <- hcl.colors(3, palette = "zissou1", alpha = 0.7)

  p1 <- density(as.numeric(mod_res$draws("indi_soc_sd")))
  p2 <- density(as.numeric(mod_res$draws("dyad_soc_sd")))
  plot(0, 0, xlim = xlim, ylim = c(0, max(p1$y, p2$y) * 1.05),
       type = "n", xaxs = "i", yaxs = "i", xlab = "estimated SD",
       ylab = "density", axes = FALSE, main = "")
  if (do_legend) {
    legend("topright",
           legend = c(labs$indi, labs$dyad, "prior"),
           col = c(xcols[1:2], "grey"), pch = c(15, 15, NA),
           lty = c(NA, NA, 1), bty = "n", cex = 0.7)
  }

  if (do_indi) polygon(p1, border = NA, col = xcols[1])
  if (do_dyadic) polygon(p2, border = NA, col = xcols[2])

  axis(1)

  pr <- brms::rstudent_t(4000, 3, 0, 1)
  pr <- density(pr, adjust = 2)
  polygon(pr, border = "grey")
  box(bty = "l")
}
