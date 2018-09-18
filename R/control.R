#' Target is a point
#'
#' @param limits
#' @param normalApprox
#' @param ctrlOverDisp
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each end
#' @param multAdj ...
#'
#' @export
pointTarget <- function(limits = 0.05, normalApprox = TRUE, crtlOverDisp = FALSE,
                        trim = 0, multAdj = "none") {
  list(method = "point", limits = limits, normalApprox = normalApprox,
       crtlOverDisp = crtlOverDisp, multAdj = multAdj, trim = trim)
}

#' Target is a distribution
#'
#' @inheritParams pointTarget
#'
#' @export
distTarget <- function(limits = 0.05, trim = 0) {
  list(method = "distribution", limits = limits, trim = trim)
}
