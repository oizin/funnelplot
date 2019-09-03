#' Example dataset for the funnelplot pacakge
#'
#' A synthetic dataset containing XX hospitals and XX patients. Hospitals are being compared on the variable
#' test, with gender and age as potential adjustment variables.
#'
#' @format A data frame with 20445 rows and 7 variables:
#' \describe{
#'   \item{hosp_id}{cluster ID}
#'   \item{patient_id}{observation ID}
#'   \item{gender}{patient gender (0:male, 1:female)}
#'   \item{age}{age of patient (standardised)}
#'   \item{test}{outcome variable: whether a procedure was successful}
#' }
"example_data"
