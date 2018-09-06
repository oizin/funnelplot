#' Hello.
#'
#' @params formula
#'
getFunnelFormula <- function(formula) {
  vars <- all.vars(formula)
  p <- length(vars)
  if(p > 2) {
    covar <- paste0(vars[2:(p-1)],collapse = "+")
    newForm <- paste(vars[1],covar,sep = "~")
  } else {
    newForm <- paste(vars[1],"1",sep = "~")
  }
  list(newForm = as.formula(newForm),
       outcome = vars[1],
       group = vars[p])
}
