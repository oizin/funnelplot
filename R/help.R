#########################################################################
# Helper functions for creating funnel plots
#
#########################################################################

#' Edit funnel plot formula
#'
#' Edit formula passed to funnel to enable call to stats::glm
#'
#' @param ff formula passed to funnel.
#'
getFunnelFormula <- function(ff) {
  # formula without cluster term
  ff <- Formula::Formula(ff)
  newff <- stats::formula(ff,lhs=1,rhs=1)
  # null formula (no covariates)
  vars <- all.vars(ff)
  p <- length(vars)
  nullff <- stats::as.formula(paste(vars[1],"1",sep = "~"))
  # return
  list(newForm = newff,
     nullForm = nullff,
     outcome = vars[1],
     id = vars[p])
}
