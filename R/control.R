#' Target is a point
#'
#' A description of how the funnel plot target and control limits should be constructed for a point target.
#'
#' @param limits a numeric vector containing the (1-limit)100\% values for the control limits.
#' @param normalApprox a logical value indicating whether use a normal approximation.
#' @param crtlOverDisp a logical value indicating whether to adjust control limits to account for overdispersion.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each end.
#' @param multAdj a character indicating which method to use to multiplicity adjustment. Default = "none", other options are
#' "fdr" and "bonferroni".
#' @param standardised a logical value indicating whether standardise the per cluster outcome. For binary data this is y = observed/expected.
#'
#' @export
pointTarget <- function(limits = 0.05, normalApprox = TRUE, crtlOverDisp = FALSE,
                        trim = 0, multAdj = "none", standardised = FALSE) {
  list(method = "point", limits = limits, normalApprox = normalApprox,
       crtlOverDisp = crtlOverDisp, multAdj = multAdj, trim = trim, standardised = standardised)
}

#' Check arguments on pointTarget
#'
#' Checks whether the arguments to pointTarget are correct.
#'
#' @param x list created from call to pointTarget.
#'
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
#' A description of how the funnel plot target and control limits should be constructed for a distribution target.
#'
#' @inheritParams pointTarget
#'
#' @export
distTarget <- function(limits = 0.05, trim = 0) {
  list(method = "distribution", limits = limits, trim = trim)
}

#' Check arguments on distTarget
#'
#' Checks whether the arguments to distTarget are correct.
#'
#' @param x list created from call to distTarget.
#'
check_distTarget <- function(x) {
  assertthat::assert_that(all(x$limits > 0.0) & all(x$limits < 0.5))
  assertthat::assert_that(x$trim >= 0.0 & x$trim < 0.5)
}
