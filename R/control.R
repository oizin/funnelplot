#' Target is a point
#'
#' @param limits ...
#' @param normalApprox ...
#' @param ctrlOverDisp ...
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each end
#' @param multAdj ...
#' @param standardised ...
#'
#' @export
pointTarget <- function(limits = 0.05, normalApprox = TRUE, crtlOverDisp = FALSE,
                        trim = 0, multAdj = "none", standardised = FALSE) {
  list(method = "point", limits = limits, normalApprox = normalApprox,
       crtlOverDisp = crtlOverDisp, multAdj = multAdj, trim = trim, standardised = standardised)
}

#' Check arguments on pointTarget
#'
#' @param x pointTarget
check_pointTarget <- function(x) {
  assertthat::assert_that(all(x$limits > 0.0) & all(x$limits < 0.5))
  assertthat::assert_that(is.logical(x$normalApprox))
  assertthat::assert_that(is.logical(x$crtlOverDisp))
  assertthat::assert_that(x$trim >= 0.0 & x$trim < 0.5)
  assertthat::assert_that(x$multAdj %in% c("none","fdr","bonferroni"))
  assertthat::assert_that(is.logical(x$standardised))
}

#' Target is a distribution
#'
#' @inheritParams pointTarget
#'
#' @export
distTarget <- function(limits = 0.05, trim = 0) {
  list(method = "distribution", limits = limits, trim = trim)
}

#' Check arguments on distTarget
#'
#' @param x distTarget
check_pointTarget <- function(x) {
  assertthat::assert_that(all(x$limits > 0.0) & all(x$limits < 0.5))
  assertthat::assert_that(x$trim >= 0.0 & x$trim < 0.5)
}
