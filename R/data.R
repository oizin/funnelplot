#' Example dataset for the funnelplot pacakge
#'
#' A synthetic dataset containing XX hospitals and XX patients. Hospitals are being compared on the variable
#' test, with X1 and X1 as potential adjustment variables.
#'
#' @format A data frame with 20445 rows and 7 variables:
#' \describe{
#'   \item{hosp_id}{cluster ID}
#'   \item{s0}{remove}
#'   \item{n}{remove}
#'   \item{patient_id}{observation ID}
#'   \item{var}{variable - rename}
#'   \item{age}{age of patient (standardised)}
#'   \item{test}{outcome variable: whether a test was performed}
#' }
"example_data"
