#########################################################################
# Helper functions for creating funnel plots
#
#########################################################################

#' Edit funnel plot formula
#'
#' Edit formula passed to funnel to enable call to stats::glm
#'
#' @param formula formula passed to funnel.
#'
getFunnelFormula <- function(formula) {
  # formula without cluster term
  formula <- Formula::Formula(formula)
  newForm <- stats::formula(formula,lhs=1,rhs=1)
  vars <- all.vars(formula)
  p <- length(vars)
  # return
  list(newForm = newForm,
     outcome = vars[1],
     id = vars[p])
}
