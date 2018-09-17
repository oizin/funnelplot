#' Hello.
#'
#' @params formula
#' @export
getFunnelFormula <- function(formula) {
  vars <- all.vars(formula)
  p <- length(vars)
  nullForm <- paste(vars[1],"1",sep = "~")
  if(p > 2) {
    covar <- paste0(vars[2:(p-1)],collapse = "+")
    newForm <- paste(vars[1],covar,sep = "~")
  } else {
    newForm <- nullForm
  }
  list(newForm = as.formula(newForm),
       nullForm = as.formula(nullForm),
       outcome = vars[1],
       id = vars[p])
}
